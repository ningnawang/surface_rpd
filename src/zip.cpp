#include "zip.hpp"

#include <zlib.h>

#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

// #####################################################################
//  Simplified version of partio
// #####################################################################
namespace ZIP {

template <class T>
inline void Write_Primitive(std::ostream& stream, const T& x) {
  stream.write(&(char&)x, sizeof(T));
}

// #####################################################################
//  class GZipFileHeader
// #####################################################################
struct GZipFileHeader {
  unsigned char magic0, magic1;  // magic should be 0x8b,0x1f
  unsigned char cm;              // compression method 0x8 is gzip
  unsigned char flags;           // flags
  unsigned int modtime;          // 4 byte modification time
  unsigned char flags2;          // secondary flags
  unsigned char os;              // operating system 0xff for unknown
  unsigned short crc16;          // crc check
  unsigned int crc32;

  GZipFileHeader()
      : magic0(0),
        magic1(0),
        flags(0),
        modtime(0),
        flags2(0),
        os(0),
        crc16(0),
        crc32(0) {}

  void Write(std::ostream& ostream) {
    magic0 = 0x1f;
    magic1 = 0x8b;
    cm = 8;
    flags = 0;
    os = 0xff;
    Write_Primitive(ostream, magic0);
    Write_Primitive(ostream, magic1);
    Write_Primitive(ostream, cm);
    Write_Primitive(ostream, flags);
    Write_Primitive(ostream, modtime);
    Write_Primitive(ostream, flags2);
    Write_Primitive(ostream, os);
  }
};

// #####################################################################
//  class ZipStreambufCompress
// #####################################################################
class ZipStreambufCompress : public std::streambuf {
  static const int buffer_size = 512;
  std::ostream& ostream;  // owned when header==0 (when not part of zip file)

  z_stream strm;
  unsigned char in[buffer_size], out[buffer_size];

  GZipFileHeader gzip_header;
  unsigned int header_offset;
  unsigned int uncompressed_size;
  unsigned int crc;

  bool valid;

 public:
  ZipStreambufCompress(std::ostream& stream) : ostream(stream), valid(true) {
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    int ret = deflateInit2(&strm, Z_BEST_SPEED, Z_DEFLATED, -MAX_WBITS, 8,
                           Z_DEFAULT_STRATEGY);
    if (ret != Z_OK) {
      std::cerr << "libz: failed to deflateInit" << std::endl;
      valid = false;
      return;
    }
    setg(0, 0, 0);
    setp((char*)in, (char*)(in + buffer_size - 4));  // we want to be 4 aligned
    // Write header
    header_offset = static_cast<unsigned int>(stream.tellp());
    gzip_header.Write(ostream);

    uncompressed_size = crc = 0;
  }

  virtual ~ZipStreambufCompress() {
    if (valid) {
      // std::cout << "Deconstruct ZipStreambufCompress" << std::endl;
      process(true);
      deflateEnd(&strm);
      Write_Primitive(ostream, crc);
      Write_Primitive(ostream, uncompressed_size);
      // use delete can guarantee compression is finished before next write &
      // prevent memory leak
      delete &ostream;
    }
  }

 protected:
  int process(bool flush) {
    if (!valid) return -1;
    strm.next_in = (Bytef*)pbase();
    strm.avail_in = static_cast<uInt>(pptr() - pbase());
    while (strm.avail_in != 0 || flush) {
      strm.avail_out = buffer_size;
      strm.next_out = (Bytef*)out;
      int ret = deflate(&strm, flush ? Z_FINISH : Z_NO_FLUSH);
      if (!(ret != Z_BUF_ERROR && ret != Z_STREAM_ERROR)) {
        valid = false;
        std::cerr << "gzip: gzip error " << strm.msg << std::endl;
        return -1;
      }
      int generated_output = static_cast<int>(strm.next_out - (Bytef*)out);
      ostream.write((char*)out, generated_output);
      if (ret == Z_STREAM_END) break;
    }
    // update counts, crc's and buffers
    int consumed_input = static_cast<int>(pptr() - pbase());
    uncompressed_size += consumed_input;
    crc = crc32(crc, (Bytef*)in, consumed_input);
    setp(pbase(), pbase() + buffer_size - 4);
    return 1;
  }

  virtual int sync() {
    if (pptr() && pptr() > pbase()) return process(false);
    return 0;
  }

  virtual int underflow() {
    std::runtime_error("Attempt to read write only ostream");
    return 0;
  }

  virtual int overflow(int c = EOF) {
    if (c != EOF) {
      *pptr() = static_cast<char>(c);
      pbump(1);
    }
    if (process(false) == EOF) return EOF;
    return c;
  }

  // assignment operator declared and not defined, to suppress warning 4512 for
  // Visual Studio
  ZipStreambufCompress& operator=(const ZipStreambufCompress& _Right);
};

// #####################################################################
//  Class ZIP_FILE_OSTREAM
// #####################################################################
//  Class needed because ostream cannot own its streambuf
class ZIP_FILE_OSTREAM : public std::ostream {
  ZipStreambufCompress buf;

 public:
  ZIP_FILE_OSTREAM(std::ostream& ostream) : std::ostream(&buf), buf(ostream) {}

  virtual ~ZIP_FILE_OSTREAM() {}
};

std::ostream* Gzip_Out(const std::string& filename, std::ios::openmode mode) {
  std::ofstream* outfile = new std::ofstream(filename.c_str(), mode);
  return new ZIP_FILE_OSTREAM(*outfile);
}
}  // namespace ZIP