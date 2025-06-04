import os
import json

import os
import json

from tqdm import  tqdm


def compare_jsons(folder1, folder2):
    # Get all json files in folder1 (subset)
    folder1_files = sorted([f for f in os.listdir(folder1) if f.endswith('.json')])
    count=0
    for file_name in tqdm(folder1_files):
        path1 = os.path.join(folder1, file_name)
        path2 = os.path.join(folder2, file_name)

        # Check if the file exists in folder 2
        if not os.path.exists(path2):
            print(f"{file_name} is missing in folder 2")
            continue

        with open(path1, 'r') as f1, open(path2, 'r') as f2:
            data1 = json.load(f1)
            data2 = json.load(f2)

            # Check if 'index' and 'pred_acgt' exist in both files
            if 'index' not in data1 or 'index' not in data2:
                print(f"Key 'index' is missing in {file_name}")
                continue
            if 'pred_acgt' not in data1 or 'pred_acgt' not in data2:
                print(f"Key 'pred_acgt' is missing in {file_name}")
                continue

            # Compare the values
            if data1['index'] != data2['index']:
                print(f"Mismatch in 'index' for {file_name}: {data1['index']} != {data2['index']}")

            if data1['pred_acgt'] != data2['pred_acgt']:
                print(f"Mismatch in 'pred_acgt' for {file_name}: {data1['pred_acgt']} != {data2['pred_acgt']}")
                label1 = data1['label_acgt']
                label2 = data2['label_acgt']
                print(f"Labels is {label1}, {label2}")
                if label1 != label2:
                    print ("Error")
                else:
                    if label1 !=data1['pred_acgt']:
                        print("new pred is wrong")
                    if label2 != data2['pred_acgt']:
                        print("old pred is wrong")
                    if label1!=data1['pred_acgt'] and label2 ==data2['pred_acgt']:
                        count =count+1
                        print("Wrong when correct")
    print("number of wrong correct")
    print(count)


def cluster_parser():
    # Path to the folder containing the original JSON files
    input_folder = '/Users/omersabary/Documents/dataNanopore-pilotSep11/inf_itai_final/dnnrobustness/full_nanopore.DNAformer_s_v2_XL/full_nanopore_size_16/'
    output_folder = '/Users/omersabary/Documents/dataNanopore-pilotSep11/inf_itai_final/dnnrobustness/full_nanopore.DNAformer_s_v2_XL/pseudo_clusters_samples/'

    # Make sure the output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Iterate over all the JSON files in the input folder
    for filename in os.listdir(input_folder):
        if filename.endswith('.json'):
            # Construct full file path
            file_path = os.path.join(input_folder, filename)

            # Open and read the JSON file
            with open(file_path, 'r') as f:
                data = json.load(f)

            # Extract the required information
            index = data['index']
            label_acgt = data['label_acgt']
            noisy_copies = data['noisy_copies']
            noisy_copies_with_index = []
            for copy in noisy_copies:
                copy_new = index + copy
                noisy_copies_with_index.append(copy_new)

            # Create the new JSON structure
            new_data = {
                "index": index,
                "data": index + label_acgt,  # Rename label_acgt to data
                "noisy_copies": noisy_copies_with_index
            }

            # Write the new JSON file, named with the value of "index"
            output_file = os.path.join(output_folder, f"{index}.json")
            with open(output_file, 'w') as outfile:
                json.dump(new_data, outfile, indent=4)

    print("All files processed successfully.")

if __name__ == "__main__":
    folder1 = "/Users/omersabary/PycharmProjects/DeepDNA-EndocderDecoder_Git/results/"  # Path to the first folder (subset)
    folder1 = "../results/"
    folder2 = "/Users/omersabary/Documents/dataNanopore-pilotSep11/inf_itai_final/dnnrobustness/full_nanopore.DNAformer_s_v2_XL/full_nanopore_size_16/"  # Path to the second folder
    folder2 = "./results_full/"  # Path to the second folder

    compare_jsons(folder1, folder2)
    #cluster_parser()

