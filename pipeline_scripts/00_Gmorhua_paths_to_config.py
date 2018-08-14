#!/usr/bin/env python

#Script specific for Gmorhua to change full paths to local_config for use with script 00.

def main():
    old_file_name = "Gmorhua_lane1_paths.txt"
    new_file_name = "Gmorhua_local_config.txt"

    old_file = open(old_file_name,"r")
    new_file = open(new_file_name,"w")

    for line in old_file:
        line = line.strip()
        if line[0] == "#":
            pass
        else:
            split_line = line.split("/")
            filename = split_line[-1]
        
            #Grab sample ID, if from sixcod30x dir, add pop info, otherwise, use existing sample name
            if split_line[7] == "sixcod30x":
                sample_num = filename.split("_")[0]
                sample_pop = filename.split("_")[1]
                sample = '%s%s'%(sample_num,sample_pop)
            else:
                sample = filename.split("_")[0]
            read_num_path = filename.split(".")[-3]
            if read_num_path == "R1":
                read_num = "1"
            elif read_num_path == "R2":
                read_num = "2"
            else:
                print("Read num not equal to R1 or R2 for %s"%line)
                
            #Create lane 2 path
            line_lane2 = line.replace("Lane1","Lane2")
        
            #Write Lane 1 line
            new_file.write('%s\t%s\t%s\t%s_Lane1\n'%(line,read_num,sample,sample))
            
            #Write Lane 2 line
            new_file.write('%s\t%s\t%s\t%s_Lane2\n'%(line_lane2,read_num,sample,sample))
    
    old_file.close()
    new_file.close()

    
if __name__ == "__main__":
    main()