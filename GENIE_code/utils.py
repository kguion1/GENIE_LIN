import pandas as pd
import numpy as np
import anndata as ad

#Download data_clinical_patient.txt, data_clinical_sample.txt, data_mutations_extended.txt
#Make sure files are in same folder
def preprocess( path_to_folder_with_files, \
                patient_file="data_clinical_patient.txt", sample_file="data_clinical_sample.txt", maf_file="data_mutations_extended.txt" ):
    """move from downloading data to necessary data frames"""

    #Read in files
    clinical_patient = pd.read_csv( "{}/{}".format(folder_with_files, patient_file), \
                                                    sep="\t", low_memory=False, skiprows=4, index_col=0  )
    clinical_sample = pd.read_csv( "{}/{}".format(folder_with_files, sample_file), \
                                                    sep="\t", low_memory=False, skiprows=4, index_col=0 )
    maf = pd.read_csv( "{}/{}".format(folder_with_files, maf_file), \
                                                    sep=",", low_memory=False, )

    


def create_anndata_mut_counts( sample_data, clinic_data, maf_data, h5ad_output_file ):
    """create anndata frame with mutation counts"""

    obs_col= list(sample_data.columns) + list(clinic_data.columns)
    obs_pandas = pd.DataFrame( index=list(sample_data.index), columns= obs_col)
    #obs_samp = pd.DataFrame( index=list(sample_data.index), list(sample_data.columns) )
    #obs_patient = pd.DataFrame( index=list(sample_data.index), list(clinic_data.columns) )
    print("Creating obs pandas...")
    for samp in sample_data.index:
        patient = sample_data.loc[samp, "PATIENT_ID"]
        row = sample_data.loc[samp, :] + clinic_data.loc[patient, :]
        obs_pandas.loc[ samp ] = row

    print("Obs pandas created!")
    print("Creating X pandas...")
    X_pandas = pd.DataFrame( index=list(sample_data.index), columns=np.unique(maf_data["Hugo_Symbol"]) )
    for g in np.unique( list(maf_data["Hugo_Symbol"]) ):
        g_df = maf_data[ maf_data["Hugo_Symbol"] == g ]
        unique_samples = np.unique( g_df["Tumor_Sample_Barcode"], return_counts=True )
        for i in range( len( unique_samples[0] ) ):
            X_pandas.loc[ unique_samples[0][i], g ] = unique_samples[1][i]

    print("X pandas created!")

    anndata = ad.AnnData( X= X_pandas, obs= obs_pandas )
    anndata.obs_names_make_unique()

    anndata.write_h5ad(h5ad_output_file)

def create_kras_mut_counts( sample_data, clinic_data, maf_data, output_file ):
    """create a pandas kras mut counts with clinical data"""

    #check that all tumor samples in maf are also in sample data
    assert all( samp in list(sample_data.index) for samp in np.unique( maf_data["Tumor_Sample_Barcode"] )  )

    obs_col= ["KRAS_COUNT"] + list(sample_data.columns) + list(clinic_data.columns)
    obs_pandas = pd.DataFrame( index=list(sample_data.index), columns= obs_col)
    #obs_samp = pd.DataFrame( index=list(sample_data.index), list(sample_data.columns) )
    #obs_patient = pd.DataFrame( index=list(sample_data.index), list(clinic_data.columns) )
    print("Creating obs pandas...")
    for samp in sample_data.index:
        patient = sample_data.loc[samp, "PATIENT_ID"]
        row = [0] + list(sample_data.loc[samp, :]) + list(clinic_data.loc[patient, :])
        assert len(row) == len(obs_pandas.columns)
        obs_pandas.loc[ samp ] = row

    print("Obs pandas created!")

    obs_pandas.sort_index(inplace=True)

    print("Add kras column...")
    #X_pandas = pd.DataFrame( index=list(sample_data.index), columns=np.unique(maf_data["Hugo_Symbol"]) )
    # for g in np.unique( list(maf_data["Hugo_Symbol"]) ):
    #     g_df = maf_data[ maf_data["Hugo_Symbol"] == g ]
    #     unique_samples = np.unique( g_df["Tumor_Sample_Barcode"], return_counts=True )
    #     for i in range( len( unique_samples[0] ) ):
    #         X_pandas.loc[ unique_samples[0][i], g ] = unique_samples[1][i]

    kras_df = maf_data[ maf_data["Hugo_Symbol"] == "KRAS" ]
    unique_samples = np.unique(kras_df["Tumor_Sample_Barcode"], return_counts=True )
    for i in range( len( unique_samples[0] ) ):
        obs_pandas.loc[ unique_samples[0][i], "KRAS_COUNT" ] = unique_samples[1][i]
    # sorted_names_indices = np.argsort( unique_samples[0] )
    # sorted_names = unique_samples[0][sorted_names_indices]
    # print(sorted_names[1:10])
    # print(list(obs_pandas.index)[1:10])
    # assert sorted_names == list(obs_pandas.index)
    # sorted_counts = unique_samples[1][ sorted_names ]
    # obs_pandas.loc[ "KRAS_Count" ] = sorted_counts


    obs_pandas.to_csv(output_file)

def create_kras_mut_counts_per_sample( kras_mut_counts_sample, clinic_data, output_file ):
    """Create kras mutation counts per patient"""

    col = ["KRAS_MUT"] + list(clinic_data.columns) + ["NUM_SAMP"] + ["SAMP_CONSIS"]
    pandas = pd.DataFrame( index=list(clinic_data.index), columns=col )

    for patient in clinic_data.index:
        kras_mut_sample_patient_p = kras_mut_counts_sample[ kras_mut_counts_sample[ "PATIENT_ID" ] == patient ]
        kras = False
        if sum( kras_mut_sample_patient_p[ "KRAS_COUNT" ] ) >= 1:
            kras = True

        samp_consis = True
        if len( np.unique(kras_mut_sample_patient_p["KRAS_COUNT"]) ) > 1:
            samp_consis = False

        num_samp = kras_mut_sample_patient_p.shape[0]

        row = [kras] + list(clinic_data.loc[patient, :]) + [num_samp] + [samp_consis]
        pandas.loc[ patient ] = row

    pandas.to_csv(output_file)


if __name__ == '__main__':
    sd = pd.read_csv("GENIE_data/GENIE_v8/genie_v8_data_clinical_sample_colorectal.csv", index_col=0)
    cd = pd.read_csv("GENIE_data/GENIE_v8/genie_v8_data_clinical_patient_colorectal.csv", index_col=0)
    maf = pd.read_csv("GENIE_data/GENIE_v8/genie_v8_data_mutations_extended_colorectal.csv", index_col=0, low_memory=False)
    #create_anndata_mut_counts( sample_data=sd, clinic_data=cd ,maf_data=maf, \
    #                           h5ad_output_file= "GENIE_data/GENIE_v8/genie_v8_anndata_mutation_count.h5ad")
    #create_kras_mut_counts( sample_data=sd, clinic_data=cd ,maf_data=maf, \
    #                             output_file= "GENIE_data/GENIE_v8/genie_v8_kras_mutation_count_sample.csv")
    kras_count = pd.read_csv("GENIE_data/GENIE_v8/genie_v8_kras_mutation_count_sample.csv", index_col=0)
    create_kras_mut_counts_per_sample( kras_mut_counts_sample = kras_count, clinic_data=cd, \
                                output_file = "GENIE_data/GENIE_v8/genie_v8_kras_mutation_count_patient.csv" )
