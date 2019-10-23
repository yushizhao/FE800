#include <iostream>
#include <algorithm>
#include <fstream>
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
#include <ctime>
#include <sstream>

#include "AnuRandom.hpp"

using namespace std;

namespace
{
// look-up table for char to hex conversion (vim rulez! ;))
const string g_lut[256] = {
                            "00", "01", "02", "03", "04", "05", "06", "07",
                            "08", "09", "0a", "0b", "0c", "0d", "0e", "0f",
                            "10", "11", "12", "13", "14", "15", "16", "17",
                            "18", "19", "1a", "1b", "1c", "1d", "1e", "1f",
                            "20", "21", "22", "23", "24", "25", "26", "27",
                            "28", "29", "2a", "2b", "2c", "2d", "2e", "2f",
                            "30", "31", "32", "33", "34", "35", "36", "37",
                            "38", "39", "3a", "3b", "3c", "3d", "3e", "3f",
                            "40", "41", "42", "43", "44", "45", "46", "47",
                            "48", "49", "4a", "4b", "4c", "4d", "4e", "4f",
                            "50", "51", "52", "53", "54", "55", "56", "57",
                            "58", "59", "5a", "5b", "5c", "5d", "5e", "5f",
                            "60", "61", "62", "63", "64", "65", "66", "67",
                            "68", "69", "6a", "6b", "6c", "6d", "6e", "6f",
                            "70", "71", "72", "73", "74", "75", "76", "77",
                            "78", "79", "7a", "7b", "7c", "7d", "7e", "7f",
                            "80", "81", "82", "83", "84", "85", "86", "87",
                            "88", "89", "8a", "8b", "8c", "8d", "8e", "8f",
                            "90", "91", "92", "93", "94", "95", "96", "97",
                            "98", "99", "9a", "9b", "9c", "9d", "9e", "9f",
                            "a0", "a1", "a2", "a3", "a4", "a5", "a6", "a7",
                            "a8", "a9", "aa", "ab", "ac", "ad", "ae", "af",
                            "b0", "b1", "b2", "b3", "b4", "b5", "b6", "b7",
                            "b8", "b9", "ba", "bb", "bc", "bd", "be", "bf",
                            "c0", "c1", "c2", "c3", "c4", "c5", "c6", "c7",
                            "c8", "c9", "ca", "cb", "cc", "cd", "ce", "cf",
                            "d0", "d1", "d2", "d3", "d4", "d5", "d6", "d7",
                            "d8", "d9", "da", "db", "dc", "dd", "de", "df",
                            "e0", "e1", "e2", "e3", "e4", "e5", "e6", "e7",
                            "e8", "e9", "ea", "eb", "ec", "ed", "ee", "ef",
                            "f0", "f1", "f2", "f3", "f4", "f5", "f6", "f7",
                            "f8", "f9", "fa", "fb", "fc", "fd", "fe", "ff"
                          };
} // unnamed namespace

int main(void)
{

  AnuRandom rnd;
  AnuRandom::Data data;
  //cout<<"got "<<data.size()<<" bytes"<<endl;
  
  //static ofstream outfile ("ANU_B", ofstream::binary | fstream::app);
  static FILE* outfile;
  outfile = fopen("ANU_25Jun2014_100MB_1+2", "ab");
  //static stringstream buffer(stringstream::out|stringstream::binary);
	static char buffer[1024]; 
	static int b = 0;
	bool success = false;
	int retry = 0;
	
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> elapsed_seconds;
	std::time_t start_time;
	
	while (!success && retry<99999999)
	{
		try{
			// web service calls
			  start = std::chrono::system_clock::now();
			  start_time = std::chrono::system_clock::to_time_t(start);
			  cout << "connecting at " << std::ctime(&start_time) << "\n";
			  
			  for (int i=0; i<1024*1024; i++){
				  //for (int j=0; j<1024; j++){
					  data=rnd.read();
					  for_each( data.begin(), data.end(), [](char c) {buffer[b]=c; b++;} );
				 // }
				//cout << i << "buffer" << "\n";
				//end = std::chrono::system_clock::now();
				//elapsed_seconds = end-start;
				//cout << "buffer elapsed time: " << elapsed_seconds.count() << "\n";
				
				//outfile.write(buffer.str().c_str(), buffer.str().length());
				//cout << sizeof(data[0]) << " " << data.size() << "\n";
				//cout << b << "\n";
				fwrite(&buffer, 1, 1024, outfile);
				b = 0;
				//end = std::chrono::system_clock::now();
				//elapsed_seconds = end-start;
				//cout << "write elapsed time: " << elapsed_seconds.count() << "\n";
				
				//buffer.str( std::string() );
				//buffer.clear();
				//outfile << '\n';
			  }
			end = std::chrono::system_clock::now();
			elapsed_seconds = end-start;
			success = true;
			cout << "success" << "\n";
			cout << "elapsed time: " << elapsed_seconds.count() << "\n";
		} catch(...) {
			end = std::chrono::system_clock::now();
			elapsed_seconds = end-start;
			cout << "elapsed time: " << elapsed_seconds.count() << "\n";
			cout << "retry... ";
			std::this_thread::sleep_for (std::chrono::seconds(60));
			retry ++;
			cout << retry << "\n";
		}
	}

  // for (int i=0; i<1024*1024; i++){
    // data=rnd.read();
    // for_each( data.begin(), data.end(), [](uint8_t c) { outfile << c;} );
    // outfile << '\n';
  // }
  // for (int i=0; i<1; i++){
  //   data=rnd.read();
  //   for_each( data.begin(), data.end(), [](uint8_t c) { outfile << g_lut[c];} );
  //   //outfile << '\n';
  // }
  //for_each( data.begin(), data.end(), [](uint8_t c) { cout<<g_lut[c]; outfile << c;} );
  //for_each( data.begin(), data.end(), [](uint8_t c) { cout<<g_lut[c]; outfile << g_lut[c];} );
  //for_each( data.begin(), data.end(), [](char c) { cout<<c; } );
  //cout<<data;
  fclose(outfile);
  //cout<<endl;
  return 0;
}
