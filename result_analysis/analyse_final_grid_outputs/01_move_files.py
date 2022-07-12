from distutils.dir_util import copy_tree
import glob
import os

SOURCE_DIR = "/beegfs/rpm/projects/monica/out"
TARGET_DIR = "/beegfs/jaenicke/klimertrag_temp/raster/sim-yields/11_final_results"

FOLDER_LIST = [
    "jaenicke_3184_2022-16-January_195527",
    "jaenicke_3183_2022-16-January_195435",
    "jaenicke_3182_2022-16-January_195402",
    "jaenicke_3181_2022-16-January_195330",
    "jaenicke_3180_2022-16-January_195300",
    "jaenicke_3179_2022-16-January_195224",
    "jaenicke_3170_2022-13-January_153227",
    "jaenicke_3169_2022-13-January_153147",
    "jaenicke_3168_2022-13-January_153115",
    "jaenicke_3167_2022-13-January_153038",
    "jaenicke_3166_2022-13-January_153009",
    "jaenicke_3165_2022-13-January_152906",
    "jaenicke_3164_2022-13-January_152831",
    "jaenicke_3163_2022-13-January_152801",
    "jaenicke_3162_2022-13-January_152728",
    "jaenicke_3161_2022-13-January_152655",
    "jaenicke_3160_2022-13-January_152616",
    "jaenicke_3159_2022-13-January_152547",
    "jaenicke_3158_2022-13-January_152445",
    "jaenicke_3157_2022-13-January_152414",
    "jaenicke_3156_2022-13-January_152343",
    "jaenicke_3155_2022-13-January_152254",
    "jaenicke_3154_2022-13-January_152206",
    "jaenicke_3153_2022-13-January_152129",
    "jaenicke_3152_2022-13-January_152048",
    "jaenicke_3151_2022-13-January_152018",
    "jaenicke_3150_2022-13-January_151942"
    ]

def move_files(folder_list, source_dir, target_dir):


    for folder in FOLDER_LIST:
        from_dirs = glob.glob(f"{source_dir}/{folder}/out/*")
        for from_dir in from_dirs:
            base_name = os.path.basename(from_dir)
            to_dir = f"{target_dir}/{base_name}"
            print(from_dir, '\n', to_dir)
            copy_tree(from_dir, to_dir)
            print("Done!")

def print_number_files_per_folder(target_dir, file_ext):
    from_dirs = glob.glob(f"{target_dir}/*")
    for from_dir in from_dirs:
        file_list = glob.glob(f"{from_dir}/*{file_ext}")
        size_list = [os.path.getsize(path) for path in file_list]
        size_list = set(size_list)
        print(from_dir, len(file_list), "\n", )


def print_crops(target_dir, file_ext):
    from_dirs = glob.glob(f"{target_dir}/*")
    crops = []        
    for from_dir in from_dirs:
        file_list = glob.glob(f"{from_dir}/*{file_ext}")
        if len(file_list) != 0:
            crop_names = [os.path.basename(file).split('_')[0] for file in file_list]
            crop_names = set(crop_names)
            crops += crop_names
    print(set(crops))

def main():
    
    move_files(folder_list=FOLDER_LIST, source_dir=SOURCE_DIR, target_dir=TARGET_DIR)

    print_number_files_per_folder(TARGET_DIR, ".asc")
    print_crops(TARGET_DIR, ".asc")

if __name__ == '__main__':
    main()