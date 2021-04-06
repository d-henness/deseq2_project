import os
import argparse
import sys

def file_is_fq(filenm):
    if filenm.endswith('.fastq') or filenm.endswith('.fastq.gz') or filenm.endswith('.fq') or filenm.endswith('.fq.gz'):
        return True
    return False

def main():
    parser = argparse.ArgumentParser("Make config file for snakemake workflows. Directories containing the string 'tmp' anywhere in their name will automatically be excluded. The user will need to manually match the directories as tumour: normal pairs")
    parser.add_argument('dir_strings', nargs = '+', help = """One or more strings that match the directories containing .fq.gz files. Example: \'python3 /path/to/make_config_2.py RG MF\' will match the directories RG7_1, RG7_PBMC, and MF5_PBMC, but not the directory RG7_1_tmp""")
    parser.add_argument('-b', '--bed_file', help = "bed_file for this run", required = True)
    parser.add_argument('-r', '--rna', action = 'store_true', help = "rna mode")
    args = parser.parse_args()

    libs = {}
    merge_libs = {}


    working_directory = os.getcwd()
    for dir_string in args.dir_strings:
        path, basename = os.path.split(dir_string)
        if path == "":
            path = working_directory
        for item in sorted(os.listdir(path)):
            if os.path.isdir(f"{path}/{item}") and item.startswith(basename):
                reads_list = []
                for filenm in sorted(os.listdir(f"{path}/{item}")):
                    if file_is_fq(filenm) and ('cut_u' not in filenm):
                        reads_list.append(f"{path}/{item}/{filenm}")
                if len(reads_list) > 0:
                    try:
                        libname = int(item)
                        libname = f"\"{libname}\""
                    except:
                        libname = f"{item}"
                    merge_libs[libname] = []
                    for i in range(int(len(reads_list) / 2)):
                        try:
                            libname_iter = int(item)
                            libname_iter = f"\"{libname_iter}_{i}\""
                        except:
                            libname_iter = f"{item}_{i}"
                        libs[libname_iter] = [reads_list[2 * i], reads_list[2 * i + 1]]
                        merge_libs[libname].append(libname_iter.strip('"'))




#    for dir_name_short in args.dir_string:
#        path, basename = os.path.split(dir_name_short)
#        if path == "":
#            path = "."
#        for i in sorted(os.listdir(path)):
#            if (os.path.isdir(f"{path}/{i}") or os.path.islink(f"{path}/{i}")) and (i.startswith(basename)) and not ('tmp' in i):
#                reads_list = []
#                try:
#                    lib = int(i)
#                    lib = f"\"{lib}\""
#                except:
#                    lib = i
#                merge_lib_string = "{:}:  [".format(lib)
#                for seq in sorted(os.listdir(f"{path}/{i}")):
#                    if (seq.endswith('.fastq') or seq.endswith('.fastq.gz') or seq.endswith('.fq.gz') or seq.endswith('.fq')) and ('cut_u' not in seq):
#                        reads_list.append(os.getcwd() + '/' + path + '/' + i + '/' + seq)
#                for j in range(int(len(reads_list)/2)):
#                    lib_string = "\"{:}_{:}\":  ['{:}', '{:}']".format(i, j, reads_list[(j * 2)], reads_list[(j * 2) + 1])
#                    libs.append(lib_string)
#                    merge_lib_string = "\"{:}{:}_{:}\", ".format(merge_lib_string, i, j)
#                merge_lib_string = "{:}]".format(merge_lib_string[:-2])
#                merge_libs.append(merge_lib_string)

    print(f"alt_bed:  {args.bed_file}")

    if args.rna:
        print("rna_libraries:")
    else:
        print("libraries:")
    for i, libname in enumerate(libs):
        print(f"  {libname}: {libs[libname]}")

    if args.rna:
        print('rna_merge_libs:')
    else:
        print('merge_libs:')
    for i, libname in enumerate(merge_libs):
        print(f"  {libname}: {merge_libs[libname]}")

if __name__ == '__main__':
    main()
