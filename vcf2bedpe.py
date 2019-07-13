#!/usr/bin/python3
import sys
#TODO fix this so it works with any vcf with any field order

if __name__ == "__main__":

    Meetup = {}
    tinder = {}
    tra_matchup = {}
    for line in sys.stdin:

        if line[0] == '#':
            continue    
        parts = line.split("\t")
        fields_arr = parts[7].split(";")
        fields = {}
        for val in fields_arr:
            fields[val.split("=")[0]] = val.split("=")[1] if len(val.split("=")) > 1 else "-"
        if parts[4] == "<DUP:TRA>":
            if fields["TRAID"] in tra_matchup:
                dl = tra_matchup[fields["TRAID"]]
                _type = "translocation" if "ISINV" not in parts[7] else "inverted-translocation"
                print( dl[0][0],
                        dl[0][1],
                        int(dl[0][1])+abs(int(dl[1]["SVLEN"])),
                        parts[0],
                        parts[1],
                        int(parts[1])+1,
                        _type,parts[9].strip(),
                        sep="\t")
                del tra_matchup[fields["TRAID"]]
            else:
                tra_matchup[fields["TRAID"]]=(parts,fields)
        elif parts[4] == "<DEL:TRA>":
            if fields["TRAID"] in tra_matchup:
                dl = tra_matchup[fields["TRAID"]]
                _type = "translocation" if "ISINV" not in dl[0][7] else "inverted-translocation"
                print( parts[0],
                        parts[1],
                        int(parts[1])+abs(int(fields["SVLEN"])),
                        dl[0][0],
                        dl[0][1],
                        int(dl[0][1])+1,
                        _type,parts[9].strip(),
                        sep="\t")

                del tra_matchup[fields["TRAID"]]
            else:
                tra_matchup[fields["TRAID"]]=(parts,fields)
        elif parts[4] == "<DUP:TANDEM>":
            print(  parts[0],
                    parts[1],
                    int(parts[1])+1,
                    parts[0],
                    int(parts[1])+int(fields["SVLEN"]),
                    int(parts[1])+1+int(fields["SVLEN"]),
                    "tandem-duplication",parts[9].strip(), 
                    sep="\t") 
        elif parts[4] == "<DUP:ISP>":
            if "END2" in fields:
                start2=fields["POS2"]
                end2=fields["END2"]
                start1=parts[1]
                end1=parts[1]+1
            else:
                start1=fields["POS2"]
                end1=fields["END"]
                start2=int(parts[1])
                end2=start2+int(fields["SVLEN"])
            _type = "duplication" if "ISINV" not in parts[7] else "inverted-duplication"
            if "CHR2" in fields:
                print(  fields["CHR2"],
                        start2,
                        end2,
                        parts[0],
                        start1,
                        end1,
                        _type,parts[9].strip(), 
                        sep="\t");
            else:
                print(  parts[0],
                        start2,
                        end2,
                        parts[0],
                        start1,
                        end1,
                        _type,parts[9].strip(), 
                        sep="\t");

        elif "SVTYPE" in fields and fields["SVTYPE"]  == "INV":
            if "SVLEN" in fields: 
                print(  parts[0],
                        parts[1],
                        int(parts[1])+1,
                        parts[0],
                        int(parts[1])+abs(int(fields["SVLEN"])),
                        int(parts[1])+1+abs(int(fields["SVLEN"])),
                        "inversion",parts[9].strip(),
                        sep="\t");
            elif "END" in fields:
                print(  parts[0],
                        parts[1],
                        int(parts[1])+1,
                        parts[0],
                        fields["END"],
                        int(fields["END"])+1,
                        "inversion",parts[9].strip(),
                        sep="\t")
        elif "SVTYPE" in fields and fields["SVTYPE"] == "DEL":
            if "SVLEN" in fields and fields["SVLEN"] != "False": 
                dels = fields["SVLEN"].split(",")
                for sdel in dels:
                    if int(sdel) < 0:
                        print(  parts[0],
                                parts[1],
                                int(parts[1])+1,
                                parts[0],
                                int(parts[1])-int(sdel),
                                int(parts[1])+1-int(sdel),
                                "deletion",parts[9].strip(),
                                sep="\t");
            elif "END" in fields:
                print(  parts[0],
                        parts[1],
                        int(parts[1])+1,
                        parts[0],
                        fields["END"],
                        int(fields["END"])+1,
                        "deletion",parts[9].strip(),
                        sep="\t")

        elif "SVTYPE2" in fields and  fields["SVTYPE2"] == "INV":
            if  fields["EVENT"] in Meetup.keys():
                ev = Meetup[fields["EVENT"]]
                if fields["MATEID"] in tinder.keys():
                    mate_pos = tinder[fields["MATEID"]]
                    print(  parts[0],
                            min(min(int(mate_pos),int(parts[1])),ev[1]),
                            max(min(int(mate_pos),int(parts[1])),ev[1]),
                            ev[0],
                            min(max(int(mate_pos),int(parts[1])),ev[2]),
                            max(max(int(mate_pos),int(parts[1])),ev[2]),
                            "inversion",
                            "rec",sep="\t")
                            
                    del tinder[fields["MATEID"]]        
                else:
                    tinder[parts[2]] = parts[1]
                del Meetup[fields["EVENT"]]
            else:
                if fields["MATEID"] in tinder.keys():
                    mate_pos = tinder[fields["MATEID"]]
                    Meetup[fields["EVENT"]] = (
                            parts[0],
                            min(int(mate_pos),int(parts[1])),
                            max(int(mate_pos),int(parts[1]))
                            )
                    del tinder[fields["MATEID"]]        
                else:
                    tinder[parts[2]] = parts[1]
                    
        else:
            pass # to implement when needed 
