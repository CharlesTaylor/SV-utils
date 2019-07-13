#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <array>
#include <cctype>
#include <vector>
#include <map>
int parse_chr( std::string &str){
	if( str[0] == 'X'){ return 22;}
	if( str[0] == 'Y'){ return 23;}
	for( auto iter = str.begin(); iter != str.end(); iter++){
		if(! isdigit(*iter)){
			std::cerr << *iter << " is not a digit in " << str << std::endl;
			return -1;
		}
	}
	return stoi(str);
}

enum class vcf_source{
	dbVar
};

enum class sv_type{
	TRA, INV, DEL, DUP, INS, OTHER,
};

struct vcf_record{
	std::string chrom;
	int pos;
	std::string id;
	std::string ref;
	std::string alt;
	std::string qual;
	std::string filter;
	std::map<std::string,std::string> info;
	vcf_record(const std::string &line){
		std::istringstream line_stream(line);
		std::string info_text;
		std::string pos_text;
		line_stream >> chrom >> pos_text >> id >> ref >> alt >> qual >> filter;
		pos = std::stoi(pos_text);
		std::string info_field;
		while(std::getline(line_stream, info_field, ';')) {
			std:size_t eq_pos = info_field.find("=");
			if( eq_pos != std::string::npos){
				std::string left = info_field.substr(0,eq_pos);
				std::string right = info_field.substr(eq_pos+1);
				info[left] = right;
			}else{
				info[info_field] = "";
			}
		}
	}
};

sv_type alt_to_type(std::string alt_string){
	if(alt_string == "<DEL>"){
		return sv_type::DEL;
	}else if(alt_string == "<DUP>"){
		return sv_type::DUP;
	}else if(alt_string == "<INV>"){
		return sv_type::INV;
	}else if(alt_string == "<INS>"){
		return sv_type::INS;
	}
	return sv_type::OTHER;
}

int dbVar(int argc, char **argv){
	for(std::string line; std::getline(std::cin,line);){
		if(line[0] == '#'){ continue; }
		vcf_record rec(line);

		sv_type type = alt_to_type(rec.alt);
		if(type == sv_type::DEL){
			int svlen = -std::stoi(rec.info["SVLEN"]);
			std::cout << rec.chrom << "\t" << rec.pos << "\t" << rec.pos+svlen << "\t" << rec.alt << "\t" << rec.info["SAMPLE"] << "\n";
		}else if(type == sv_type::DUP){
			int svlen = std::stoi(rec.info["SVLEN"]);
			std::cout << rec.chrom << "\t" << rec.pos << "\t" << rec.pos+svlen << "\t" << rec.alt  << "\t" << rec.info["SAMPLE"] << "\n";
		}else if(type == sv_type::INV){
			std::cout << rec.chrom << "\t" << rec.pos << "\t" << rec.info["END"] << "\t" << rec.alt << "\t" << rec.info["SAMPLE"] << "\n";
		}

	}	

	return 0;
}

int main(int argc, char **argv){
	if( argc == 1){
		std::cerr << "dbVar" << std::endl;
		return -1;
	}
	std::string tol(argv[1]);
	if(tol == "dbVar"){
		return dbVar(argc-1,argv+1);
	}
	std::cerr << "dbVar" << std::endl;

	return -1;
}
