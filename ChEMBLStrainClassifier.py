#!/usr/bin/env python3
import bacdive
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from tqdm import tqdm

# put your details here after creating an account with BacDive
yourusername = "user_name"
password = "your_password"
client = bacdive.BacdiveClient(yourusername,password) 



def handle_searchterms(num_dict):
    """
    Args:
        A dictionary of the format:
        (full_repository_ascession, repository_number) : BacDive_ID

    Returns:
        A dictionary of the format:
        search_term : BacDive_ID

    -- Implementation -- 
    A number of assumptions are outlined below. 

    IF a repository_number is LESS than 4 characters long
        searching by the number alone risks it matching with unrelated strains 
        e.g., BR 13 strain could potentially match names like Isolate 13
        In this case we search using the preceeding string aswell (e.g., BR 13) 
    
    IF a repository_number already appears 
    IT is ASSUMED that numbers longer than 4 characters are unlikely to be used multiple times, however in the case that th
    """
    # Initialise a dictionary of patterns and identifiers 
    pattern_dict = {}

    # create a list of numbers that repeated within the dictionary. 
    numbers_only = [item[1] for item in num_dict.keys()]
    counter = Counter(numbers_only)
    repeated_numbers = [item for item,count in counter.items() if count>1]

    # iterate through each number in the dataframe 
    for name_number in num_dict.keys():
        # 1.  IF the number is less than 4 digits, its probably best to search using its full name instead
        if len(name_number[1]) < 4:
            # Sometimes the strain names are given with no space between the asscession code and the number e.g., BR13
            # or a hyphen/colon ATCC-700603, ATCC:700603
            # so if possible, try to include an optional space, hyphen or colon

            try:
                pattern = [re.escape(item) for item in name_number[0].split()]  # 1. split up the pattern, and run escape on it to prevent errors. 
                pattern = r"[\s\-:]?".join(pattern)+r"(?!\d)"                   # 2. specify that the numbers-name bit might be interspersed by a space , a hyphen - or a colon :
                                                                                # also use negative lookahead to make sure that we dont match substrings of larger numbers
                re.compile(pattern)   # <- checks that this pattern is valid
                  
            except:
                # This didnt raise an exception any of the times i tried it, but just in case, tell it to just search using the normal name given. 
                pattern = re.escape(name_number[0])
                print(f"{pattern} could not be given optional spaces in its pattern")

            pattern_dict[pattern] = num_dict[name_number]    # <-- appends to the dictionary pattern:ascession code. 

        # 2. IF the same number appears more than once in the dictionary, it is best to search using its full name. 
        elif name_number[1] in repeated_numbers:
            try:
                pattern = [re.escape(item) for item in name_number[0].split()]  # 1. split up the pattern, and run escape on it to prevent errors. 
                pattern = r"[\s\-:]?".join(pattern)+r"(?!\d)"                   # 2. specify that the numbers-name bit might be interspersed by a space , a hyphen - or a colon :
                                                                                # also use negative lookahead to make sure that we dont match substrings of larger numbers
                re.compile(pattern)   # <- checks that this pattern is valid
                 
            except:
                # This didnt raise an exception any of the times i tried it, but just in case, tell it to just search using the normal name given. 
                pattern = re.escape(name_number[0])
                print(f"{pattern} could not be given optional spaces in its pattern")
            
            pattern_dict[pattern] = num_dict[name_number]    # <-- appends to the dictionary pattern:ascession code. 

        
        # 3. If none of the above are satisfied, it is assumed that the number is safe to use...
        else:
            pattern = re.escape(name_number[1])         
            pattern =  r"(?<!\d)"+pattern+r"(?!\d)"   # <- again prevent substrings from matching. 
            pattern_dict[pattern] = num_dict[name_number]
    return pattern_dict


def culture_col_lookup(searchterm):
    """
    Given a species name queries baccdive and returns a lookup dictionary
    which maps the culture collection name to the BacDive ID of a species
    allowing us to group the species names. 
    """

    culture_dive = {} 
    client.search(taxonomy=searchterm)

    # .retrieve gives us a list of nested dictionaries from which we extract the useful data   
    for strain in client.retrieve():

        # Extract bacdive ID for that strain
        bacdive_id = strain["General"]["BacDive-ID"]

        # Extract the list of culture collection numbers associated with that strain. 
        # For a few IDs, there was no collection ID. this is not useful for our strain classifier function so these are discarded. 
        culture_collections = strain["External links"].get("culture collection no.", None)
        if culture_collections:
            # the list of names is given as a string delimeted by , so we split it to get an actual list.
            cul_col_namelist =  culture_collections.split(",")
            cul_col_namelist = [name.strip() for name in cul_col_namelist]  # also remove whitespace

            # Append all the culture collection names as keys allowing us to look up their baccdive ID quickly.
            for name in cul_col_namelist:
                culture_dive[name] = bacdive_id

    return culture_dive

def standardise_strains(dataframe, straincolumn = "assay_strain", species="Klebsiella pneumoniae"):
    """
    -- Explanation -- 
    Returns a dictionary of search patterns to group bacterial strains in the ChEMBL database. 

    Issues:
        Names input by users, when present typically refer to strains by the culture collection/repository they were acquired from 
        However even these names can be rather inconsistent.  
            e.g., ATCC 700603, 700603, ATCC700603
            To overcome this, my function preferentially searches for strains by ascession codes where possible, 
            but will include repository names where necessary to prevent duplicate values. 

        Furthermore, different repositories may house the same strain, but under different names 
            e.g., ATCC 700603 and LMG 20218 are both derived from the same source.
            To overcome this where possible, strains which are in the BacDive database are grouped by their BacDive ID to increase grouping of the data.
            Unfortunately, BacDive is not comprehensive, so this isnt possible for all strains.   


    -- Implementation-- 
    0. Preprocessing
        - Print statement informs the user how much of the dataframe simply contains NaN values. 
        - remove names which contain ";" as this suggests multiple strains (e.g., klebsiella there was one item named "NCIM 2957; ATCC 11298". 
          In these cases it is unclear which strain the MIC value came from)

    1. Using data from the BacDive database, extract the repository numbers.
        - Using the sort_numbers function we remove numbers that would be potentially problematic 
        (i.e., 13 from the BR 13 strain could potentially match names like Isolate 13, 
        so in this case it is safer to search using the full name: "BR 13")
        - Once we have filtered and created the appropriate search expressions then we use these to group names. 
        

    2. Identify strains that are associated with strain repositories (see list of repository abbreviations below)
        - Extract their catalogue numbers/identifiers i.e., '11298' from 'ATCC 11298'
        - Use these catalogue numbers to search the rest of the values. 
        This approach is used because many rows contain just the catalogue number
        or the catalogue number and the repository fused (e.g.,ATCC 700603, 700603, ATCC700603 which all clearly refer to the same strain)

    3. Identify Clinical Isolates
        - search for the term 'isolate' (not case sensitive)

    4. (for later... add an option to give custom search terms which arent automatically searched for... for example a widely used line BW25113 from Barry L. Wanner's lab wont be identified by my current methods, however may be important)
    """

    # 0. Preprocessing: 
    # create a list of unique names in the ``strain`` column
    rows_with_strains = dataframe.dropna(subset=straincolumn)
    unique_names = list(set(rows_with_strains[straincolumn].astype(str)))

    # Inform the user how many strains cant be identified. 
    print(f"Warning: {round(((len(dataframe)-len(rows_with_strains))/len(dataframe))*100, 1) }% of values have no strain entry")

    # ---- 1. BacDive strain names ---- # 
    # Creates a lookup dict which maps culture collection names to a bacdive ID
    # allowing us to group strains kept by multiple repositories
    culture_dive = culture_col_lookup(species)

    # -- Identify the ascession numbers from these culture collection names. 
    
    # <explanation of regex capturing pattern>
    #       d+ specifies must start with 1 or more digits
    #       (\.\d)? specifies it is an optional group (I later added - as some of the numbers also have dashes in them....)
    #       we capture the outermost capturing group (group 1)
    pattern = r"(\d+(\.\d+|-\d+)?[a-zA-Z]?)"    
    numbers_dict = {(collection_name, re.search(pattern,collection_name).group(1)): bacdive_id 
               for collection_name, bacdive_id in culture_dive.items() 
               if bool(re.search(pattern,collection_name)) }
    
    # -- With the numbers extracted, process these into usable patterns which minimise the chances of matching nonspecifically. 
    pattern_dict1 = handle_searchterms(numbers_dict)

    # Now we have made search patterns we can use them to build a lookup dictionary mapping 
    # the names to a BacDive number. we do this by checking each pattern against 
    # each name in the database to see if we get a match. 
    strain_lookup = {}
    for pattern in pattern_dict1.keys():
        for name in unique_names:
            if re.search(pattern,name):
                # If a pattern matches a name.....            
                if name in strain_lookup:  # if that name was already matched - thats a problem so let the user know.
                    print(f"Warning <1>: {pattern} matched {name} which was previously identified as {strain_lookup[name]}")

                strain_lookup[name] = pattern_dict1[pattern]
    
    # 2. Identify strains which are likely isolates
    # This is assuming that if an inputted name has the word isolate in it, 
    # it was probably taken deliberately for that study, and is not kept by any repository
    for name in unique_names:
        if bool(re.search("isolate", name, flags=re.IGNORECASE)):
            if name in strain_lookup:
                print(f"Warning <2>: {name} was previously identified as : {strain_lookup[name]}")

            # if Isolate in name, always name it as an isolate anyway     
            strain_lookup[name] = "Isolate"

    # 3. Regex to extract additional strain numbers
    # If an identifier for that name is not in bacDive 
    # and the name does not have the word "Isolate" then we try to extract 
    # Asscession numbers to group strains by in this subset.
    recognised_names = list(strain_lookup.keys())
    unrecognised_names = [name for name in unique_names if name not in recognised_names]


    # Using a list of repository names.  
    # Check to see which repositories are present, to save time with searching. 
    repo_list = ['IMET', 'IFO', 'CGMCC', 'CCRC', 'NCIM', 'ATCC', 'NCDO', 'LMG', 'KCTC', 'CCM', 'NBRC', 'IMG', 'CNCTC', 'DSM', 'NRRL', 'NCTC', 'CCUG', 'IAM', 'CECT', 'BCRC', 'WDCM', 'JCM', 'VTT', 'CDC', 'ECOR', 'NCDC', 'BGSC','MTCC']
    strain_longstring = " ".join(unrecognised_names)
    repos_present = [repo for repo in repo_list if bool(re.search(repo, strain_longstring))]

    numbers_dict = {}
    for repo in repos_present:
        # (I put a \b at the start to prevent NCDC from matching the same ones as CDC)
        pattern = fr"\b{repo}.*?(\d+(\.\d+|-\d+)?)"

        # creating a similar numbers dict to earlier so it can be passed into the handle_searchterms function for cleaning
        # what we need to do is to create a dictionary 
        for name in unrecognised_names:
            if bool(re.search(pattern,name)):
                # extract the ascession number
                number = re.search(pattern,name).group(1)

                # repository_name + ascession number
                fullname = f"{repo} {number}"
                # bundle this into a format that allows it to be passed into handle_searchterms 
                # and append to numbers_dict which is then used to form another pattern_dict2
                numbers_dict[(fullname,number)] = fullname


    pattern_dict2 = handle_searchterms(numbers_dict)

    for pattern in pattern_dict2.keys():
            for name in unrecognised_names:
                if re.search(pattern,name):
                    # If a pattern matches a name.....            
                    if name in strain_lookup:  # if that name was already matched - thats a problem so let the user know.
                        print(f"Warning <3>: {pattern} matched {name} which was previously identified as {strain_lookup[name]}")

                    strain_lookup[name] = pattern_dict2[pattern]
    
    # 4. Allow user to add their custom regex patterns to recognise additional strains which may not be recognised... 
    # <insert here?>

    # 5. Create the new column in the dataframe.
    dataframe["func_assigned_strain"] = dataframe[straincolumn].apply(lambda x: strain_lookup.get(x,"Unclassified"))

    #return pattern_dict1,pattern_dict2,culture_dive
    # 6. use the descriptions column to identify strains which could not be identified from the "assay_strain column"
    # Filter for items classed as unclassified so that we can now search their descriptions.
    unidentified_rows =  dataframe[dataframe["func_assigned_strain"]=="Unclassified"]

    # progress bar used on this bit as it can take a while...
    for index, row in tqdm(unidentified_rows.iterrows(), total=unidentified_rows.shape[0]):
        description = row["description"]
        if re.search("isolate", description, flags=re.IGNORECASE):
            dataframe.at[index,"func_assigned_strain"] = "Isolate"

        else:
            for pattern in pattern_dict1.keys():
                if re.search(pattern,description):
                    # check if something has already been assigned for that row.. and warn the user:
                    if dataframe.at[index,"func_assigned_strain"] != "Unclassified":
                        print(f"Warning <6> : {description} was identified by {pattern} but was previously identified by another pattern")
                        
                    dataframe.at[index,"func_assigned_strain"] = pattern_dict1[pattern]

    return dataframe





def main():
    kleb_data = pd.read_csv("C:\\Users\\Tim\\OneDrive\\__BHLS\\SSRA\\SSRA\\Data&Files\\Datasets\\Kleb_data_r.csv")
    kleb_strain = standardise_strains(dataframe=kleb_data)
    kleb_strain[["func_assigned_strain","assay_strain","description"]].drop_duplicates().to_csv("Examine_results.csv")


if __name__ == "__main__":
    main()


