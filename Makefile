mmr: mmr.cpp config.h config.cpp Alignment.h Alignment.cpp BatchData.cpp  BatchData.h  GeneralData.h OnlineData.cpp  OnlineData.h  Segment.h  Segments.cpp  Segments.h  Utils.cpp  Utils.h
	g++ -Wall -O9 -D_THREAD_SAFE -g -std=c++0x -pthread -o mmr mmr.cpp config.cpp Alignment.cpp BatchData.cpp OnlineData.cpp Segments.cpp Utils.cpp
#	g++ -Wall -g -std=c++0x -pthread -o mmr mmr.cpp config.cpp Alignment.cpp BatchData.cpp OnlineData.cpp Segments.cpp Utils.cpp

mmrX: mmr.cpp config.h config.cpp Alignment.h Alignment.cpp BatchData.cpp  BatchData.h  GeneralData.h OnlineData.cpp  OnlineData.h  Segment.h  Segments.cpp  Segments.h  Utils.cpp  Utils.h
	g++ -Wall -O9 -g -std=c++0x -pthread -o mmrX mmr.cpp config.cpp Alignment.cpp BatchData.cpp OnlineData.cpp Segments.cpp Utils.cpp
