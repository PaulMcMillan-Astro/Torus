
#include <iostream>
#include <vector>
#include "PJMebf.h"


int Ebf_GetTagNames(const std::string& filename1, 
		    std::vector<std::string> &tagnames)
{
  using namespace std;
  ebf::ebfc::EbfHeader ebfh;
  tagnames.resize(0);
  FILE* fd = fopen(filename1.c_str(), "rb");
  if (fd != NULL)
    {
      fseek(fd, 0, SEEK_END);
      int64_t loc = ftell(fd);
      fseek(fd, 0, SEEK_SET);
      while (ftell(fd) < loc)
	{
	  ebf::ebfc::EbfHeader_read(&ebfh,fd);
	  fseek(fd,ebfh.capacity_,SEEK_CUR);
	  tagnames.push_back(ebfh.name);
	  if(ebfh.ecode!=0)
	    {
	      cout<<"Ebf error: in reading header"<<endl;
	      fclose(fd);
	      return 1;
	    }
	}
      fclose(fd);
    }
  else
    {
      cout<<"Ebf error: file not found"<<endl;
    }
  return 0;
}


int Ebf_GetToriNames(const std::string& filename, 
		      std::vector<std::string> &tori) {
  
  std::vector<std::string> stuff; // all field names
  ebf::EbfDataInfo dinfo;
  if(Ebf_GetTagNames(filename, stuff)) return 1;
  for(int i=0;i!=stuff.size();i++) {
    if(stuff[i].size()>8)
      if(stuff[i].substr(stuff[i].size()-7,7)=="actions") 
	tori.push_back(stuff[i].substr(1,stuff[i].size()-9));
  }
  if(tori.size()) return 0;
  return 1;
}
