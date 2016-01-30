# This file is part of MMR, the Read Multi-Mapper Resolution tool.
#
#  MMR is free software: you can redistribute it and/or modify it
#  under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This software is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  A copy of the GNU General Public License is distributed with 
#  MMR (file LICENSE).  If not, see <http://www.gnu.org/licenses/>.
#
#  Written 2010-2015 by 
#    Andre Kahles <akahles@cbio.mskcc.org>
#    Jonas Behr <jonas_behr@web.de>
#    Gunnar R\"atsch <raetsch@cbio.mskcc.org>
# 
#  This work was funded by the Max Planck Society,
#  the German Research Foundation (DFG RA1894/2-1)
#  and Memorial Sloan Kettering Cancer Center.
#

mmr: mmr.cpp config.h config.cpp Alignment.h Alignment.cpp BatchData.cpp  BatchData.h  GeneralData.h OnlineData.cpp  OnlineData.h  Segment.h  Segments.cpp  Segments.h  Utils.cpp  Utils.h
	g++ -Wall -O9 -D_THREAD_SAFE -g -std=c++0x -pthread -o mmr mmr.cpp config.cpp Alignment.cpp BatchData.cpp OnlineData.cpp Segments.cpp Utils.cpp

debug: mmr.cpp config.h config.cpp Alignment.h Alignment.cpp BatchData.cpp  BatchData.h  GeneralData.h OnlineData.cpp  OnlineData.h  Segment.h  Segments.cpp  Segments.h  Utils.cpp  Utils.h
	g++ -Wall -g -std=c++0x -pthread -o mmr_dbg mmr.cpp config.cpp Alignment.cpp BatchData.cpp OnlineData.cpp Segments.cpp Utils.cpp

test:
	cd examples; ./test_mmr_small.sh

bigtest: 
	cd examples; ./test_mmr_large.sh
