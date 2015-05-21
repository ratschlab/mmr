/*

  This file is part of MMR, the Read Multi-Mapper Resolution tool.

  MMR is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This software is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  A copy of the GNU General Public License is distributed with 
  MMR (file LICENSE).  If not, see <http://www.gnu.org/licenses/>.

  Written 2010-2015 by 
    Andre Kahles <akahles@cbio.mskcc.org>
    Jonas Behr <jonas_behr@web.de>
    Gunnar R\"atsch <raetsch@cbio.mskcc.org>
 
  This work was funded by the Max Planck Society,
  the German Research Foundation (DFG RA1894/2-1)
  and Memorial Sloan Kettering Cancer Center.

*/

#ifndef __SEGMENT_H__
#define __SEGMENT_H__

using namespace std;

class Segment {

    public:

    unsigned long start;
    unsigned int length;
    unsigned char chr;
    unsigned char strand;
    float expectation;

    Segment(const unsigned long a, const unsigned int b, const unsigned char c,  const unsigned char d, const float e) :
        start(a), length(b), chr(c), strand(d), expectation(e) {}

    Segment():
        start(0), length(0), chr(1), strand('+'), expectation(0.0) {}

    ~Segment() {};
    
    bool operator==(const Segment &other) {
        return (start == other.start && length == other.length && chr == other.chr && strand == other.strand);
    };

};
#endif
