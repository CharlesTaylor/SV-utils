#include <iostream>
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

struct interval{

	int chr;
	int start;
	int end;
	interval (){}
	
	interval( int chr, int start, int end) :
		 chr(chr),
		 start(start),
		 end(end)
	{}
interval( std::string chr, std::string start, std::string end) :
		 chr(parse_chr(chr)),
		 start(std::stoi(start)),
		 end(std::stoi(end))
	{}

	friend std::ostream &operator<<(std::ostream &out, const interval &c){
		std::string chr = std::to_string(c.chr);
		if( c.chr == 23) { chr = "X";}
		else if( c.chr == 24){ chr = "Y";}
		out << chr << "\t" << c.start << "\t" << c.end;
		return out;
	}

	int inner_distance( const interval &other) const {
		if( start < other.start){
			return other.start - end;
		}
		return start - other.end;
	}
};

enum class tool{
	delly, lr, naibr, grocsvs,
};

enum class type{
	FF,FR,RF,RR
};

enum class sv{
	TRA, INV, DEL,
};


type stotype(std::string s1, std::string s2){
	if( s1 == "+"){
		if(s2 == "+"){
			return type::FF;
		}else if(s2 == "-"){
			return type::FR;
		}
	}else if(s1 == "-"){
		if(s2 == "+"){
			return type::RF;
		}else if(s2 == "-"){
			return type::RR;
		}
	}
	std::cerr << "Unknown RP_TYPE : " << s1 << s2 << std::endl;
	exit(-1);
}
type stotype(std::string &type_str){
	if( type_str == "TRANS_RR"){
		return type::RR;
	}else if( type_str == "TRANS_FR"){
		return type::FR;
	}else if( type_str == "TRANS_RF"){
		return type::RF;
	}else if( type_str == "TRANS_FF"){
		return type::FF;
	}else {
		std::cerr << "Unknown RP_TYPE (single) : " << type_str << std::endl;
		exit(-1);
	} 
}

type type_complement( type t){
	switch(t){
	case type::FF:
		return type::RR;
	case type::RF:
		return type::FR;
	case type::FR:
		return type::RF;
	case type::RR:
		return type::FF;
	default:
		std::cerr << "Unknown type of ordinal: " << static_cast<int>(t) << std::endl;
		exit(-1);
	}
}


std::vector<std::string> str_split(std::string &str, std::string delim){
    std::vector<std::string> splits;
    size_t p1 = 0;
    size_t p2 = 0;
    while((p2= str.find(delim,p1)) != std::string::npos){
        splits.push_back(str.substr(p1,p2-p1));
        p1 = p2+1;
    }
    splits.push_back(str.substr(p1));
    return splits;
}

//22      25915850        25915856        22      50913958        50913961        call_946        63      .       .       .       ALLELIC_FRAC=0.185779816514;BLACK1=.;BLACK2=.;BLACK_DIST1=inf;BLACK_DIST2=inf;BLACK_FRAC=0.0;FRAC_HAP_SUPPORT=0.777777777778;HAPS=0,0;HAP_ALLELIC_FRAC=1.0;LR=507.295437289;MATCHES=.;NPAIRS=11;NSPLIT=0;ORIENT=-+;PS1=23830872;PS2=50390144;RP_LR=135;RP_TYPE=TRANS_FR;SEG_DUP=.;SOURCE=SV;TYPE=UNK;ZS=HET
class bedpel {

public:
	interval src;
	interval tgt;
	type t;
	std::map<std::string,std::string> info;
	int matcc = 0;
	bedpel(tool tol, std::string &line){

		std::vector <std::string> fields;

		switch(tol){
		case tool::lr:
			fields = str_split(line, "\t");
			src = interval(fields[0],fields[1],fields[2]); 
			tgt = interval(fields[3],fields[4],fields[5]); 	
			for(std::string field :str_split(fields[11],";")){
				std::vector <std::string> pair = str_split(field,"=");
				info[pair[0]] = pair[1];
			}
	//		std::cout << line << std::endl;
			t = stotype( info["RP_TYPE"]);
			break;
		case tool::delly:
			fields = str_split(line, "\t");
			t = stotype( fields[8], fields[9]);
			src = interval(fields[0],fields[1],fields[2]); 
			tgt = interval(fields[3],fields[4],fields[5]); 	
			break;
		case tool::naibr:
			fields = str_split(line, "\t");
			src = interval(fields[0],fields[1],fields[1]);
			tgt = interval(fields[2],fields[3],fields[3]);

			t = stotype(fields[6].substr(0,1),fields[6].substr(1,1));
			break;
		case tool::grocsvs:
			fields = str_split(line, "\t");
			src = interval(fields[0],fields[1],fields[1]);
			tgt = interval(fields[2],fields[3],fields[3]);

			t = stotype(fields[4].substr(0,1),fields[4].substr(1,1));

			break;
		}
	}

	friend std::ostream &operator<<(std::ostream &out, const bedpel &c){
		out << c.src << "\t" << c.tgt << "\t" <<  static_cast<int>(c.t);
		return out;
	}
	template <sv SV>
	int agree(bedpel &other){
		switch( SV){
		case sv::TRA:
			if( src.chr == tgt.chr){return 0;}
            if(other.src.chr != src.chr){return 0;}
            if(other.tgt.chr != tgt.chr){return 0;}
            if(t != type_complement(other.t)){return 0;}
			if( src.inner_distance(other.src) < 75000){
				matcc++;
				other.matcc++;
				return 1;
			}
			if( tgt.inner_distance(other.tgt) < 75000){
				matcc++;
				
				other.matcc++;
				return 2;
			}
			return 0;
		case sv::INV:
			if(src.chr != other.src.chr){
				return 0;
			}
			if(tgt.chr != other.tgt.chr){
				return 0;
			}
			if(!(t ==  type::FF and other.t == type::RR))
			{return 0;}
			if( src.inner_distance(other.src) < 50000 &&
				tgt.inner_distance(other.tgt) < 50000){
			
				other.matcc++;
				matcc++;
				return 1;
			}
		}
        return 0;
	}
};


class call{
	interval src;
	interval tgt;
	int isinv;
	sv t;
public:
	call(bedpel b1, bedpel b2, int type, sv s) : isinv(0), t(s){
		if(b1.t == type::RR || b1.t == type::FF){
			isinv = 1;
		}
		if(type == 1){
			tgt = interval(b1.src.chr,b1.src.end,b2.src.start);
			src = interval(b1.tgt.chr,b1.tgt.start,b2.tgt.end);
		}else{
			tgt = interval(b1.tgt.chr,b1.tgt.end,b2.tgt.start);
			src = interval(b1.src.chr,b1.src.start,b2.src.end);
		}
	}
	friend std::ostream &operator<<(std::ostream &out, const call &c){
		if( c.t==sv::TRA){
			out << c.src << "\t" << c.tgt << "\t" <<  (c.isinv?"inverted-translocation":"translocation");
		}
		else if (c.t == sv::INV){
			out << c.src << "\t" << c.tgt << "\t" <<  ("inversion");
		}
		return out;
	}
};
#define CHR_CNT 24 
int grocsvs(int argc, char **argv){
	//Template hell
	std::array< std::array< std::vector<bedpel>, CHR_CNT>, CHR_CNT> bedpes;

	for(std::string line; std::getline(std::cin,line);){
		if(line[0] == '#'){ continue; }
		bedpel bp(tool::grocsvs,line);
		if(bp.src.chr ==  -1 || bp.tgt.chr == -1 ){ continue;}
		bedpes[bp.src.chr][bp.tgt.chr].push_back(bp);
	} 
	std::vector<call> matches;
	//Bracket hell
	for(int i = 0; i < bedpes.size(); i++){
		for(int j =0; j < bedpes[i].size(); j++){
			for(auto& ij : bedpes[i][j]){
				//if(ij.matcc !=0){
				//	break;
				//}
				for(auto& ji : bedpes[i][j]){
					//if(ji.matcc !=0){
					//	break;
					//}
					if(i == j && &ij == &ji){ continue;}
					int type;
				
					if( type = ij.agree<sv::INV>(ji)){						
						matches.push_back(call(ij,ji,type,sv::INV));
					}
					if( type = ij.agree<sv::TRA>(ji)){						
						matches.push_back(call(ij,ji,type,sv::TRA));
					}
				}
			}
		}
	}

	for(int i = 0; i < bedpes.size(); i++){
		for(int j =0; j < bedpes[i].size(); j++){
			for(const auto& ij : bedpes[i][j]){
				std::cerr << ij.matcc << "\t" << ij << std::endl;
			}
		}
	}
	for(auto iter = matches.begin(); iter != matches.end(); iter++){
		std::cout << *iter << std::endl;
	}
	return 0;
}

int naibr(int argc, char **argv){
	//Template hell
	std::array< std::array< std::vector<bedpel>, CHR_CNT>, CHR_CNT> bedpes;

	for(std::string line; std::getline(std::cin,line);){
		if(line[0] == '#'){ continue; }
		bedpel bp(tool::naibr,line);
		if(bp.src.chr ==  -1 || bp.tgt.chr == -1 ){ continue;}
		bedpes[bp.src.chr][bp.tgt.chr].push_back(bp);
	} 
	std::vector<call> matches;
	//Bracket hell
	for(int i = 0; i < bedpes.size(); i++){
		for(int j =0; j < bedpes[i].size(); j++){
			for(auto& ij : bedpes[i][j]){
				if(ij.matcc !=0){
					break;
				}
				for(auto& ji : bedpes[i][j]){
					if(ji.matcc !=0){
						break;
					}
					if(i == j && &ij == &ji){ continue;}
					int type;
				
					if( type = ij.agree<sv::INV>(ji)){						
						matches.push_back(call(ij,ji,type,sv::INV));
					}
					if( type = ij.agree<sv::TRA>(ji)){						
						matches.push_back(call(ij,ji,type,sv::TRA));
					}
				}
			}
		}
	}

	for(int i = 0; i < bedpes.size(); i++){
		for(int j =0; j < bedpes[i].size(); j++){
			for(const auto& ij : bedpes[i][j]){
				std::cerr << ij.matcc << "\t" << ij << std::endl;
			}
		}
	}
	for(auto iter = matches.begin(); iter != matches.end(); iter++){
		std::cout << *iter << std::endl;
	}
	return 0;
}
int longranger(int argc, char **argv){
	//Template hell
	std::array< std::array< std::vector<bedpel>, CHR_CNT>, CHR_CNT> bedpes;

	for(std::string line; std::getline(std::cin,line);){
		if(line[0] == '#'){ continue; }
		bedpel bp(tool::lr,line);
		if(bp.src.chr ==  -1 || bp.tgt.chr == -1 ){ continue;}
		bedpes[bp.src.chr][bp.tgt.chr].push_back(bp);
	} 
	std::vector<call> matches;
	//Bracket hell
	for(int i = 0; i < bedpes.size(); i++){
		for(int j =0; j < bedpes[i].size(); j++){
			for(auto& ij : bedpes[i][j]){
				for(auto& ji : bedpes[i][j]){
					if(i == j && &ij == &ji){ continue;}
					int type;
					if( type = ij.agree<sv::TRA>(ji)){
						matches.push_back(call(ij,ji,type,sv::TRA));
					}
				}
			}
		}
	}

	for(int i = 0; i < bedpes.size(); i++){
		for(int j =0; j < bedpes[i].size(); j++){
			for(const auto& ij : bedpes[i][j]){
				std::cerr << ij.matcc << "\t" << ij << std::endl;
			}
		}
	}
	for(auto iter = matches.begin(); iter != matches.end(); iter++){
		std::cout << *iter << std::endl;
	}
	return 0;
}

int delly_vcf2bedpe(int argc, char **argv){

	std::array< std::array< std::vector<bedpel>, CHR_CNT>, CHR_CNT> bedpes;

	for(std::string line; std::getline(std::cin,line);){
		if(line[0] == '#'){ continue; }
		bedpel bp(tool::delly,line);
		if(bp.src.chr ==  -1 || bp.tgt.chr == -1 ){ continue;}
		bedpes[bp.src.chr][bp.tgt.chr].push_back(bp);
	} 
	std::vector<call> matches;
	//Bracket hell
	for(int i = 0; i < bedpes.size(); i++){
		for(int j =0; j < bedpes[i].size(); j++){
			for(auto& ij : bedpes[i][j]){
				for(auto& ji : bedpes[i][j]){
					if(i == j && &ij == &ji){ continue;}
					int type;
					if( type = ij.agree<sv::TRA>(ji)){
						matches.push_back(call(ij,ji,type,sv::TRA));
					}
				}
			}
		}
	}

	for(int i = 0; i < bedpes.size(); i++){
		for(int j =0; j < bedpes[i].size(); j++){
			for(const auto& ij : bedpes[i][j]){
				std::cerr << ij.matcc << "\t" << ij << std::endl;
			}
		}
	}
	for(auto iter = matches.begin(); iter != matches.end(); iter++){
		std::cout << *iter << std::endl;
	}
	return 0;
}

int main(int argc, char **argv){
	if( argc == 1){
		std::cerr << "longranger, delly, naibr or grocsvs" << std::endl;
		return -1;
	}
	std::string tol(argv[1]);
	if(tol == "longranger"){
		return longranger(argc-1,argv+1);
	}else if (tol == "delly"){
		return delly_vcf2bedpe(argc-1,argv+1);
	
	}else if (tol =="naibr"){
		return naibr(argc-1, argv+1);
	}
	else if (tol =="grocsvs"){
		return grocsvs(argc-1, argv+1);
	}

	std::cerr << "longranger, delly, naibr or grocsvs" << std::endl;

	return -1;
}
