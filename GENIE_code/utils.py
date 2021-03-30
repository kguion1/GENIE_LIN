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
    """Create kras mutation counts per sample"""

    col = ["KRAS_MUT"] + list(clinic_data.columns) + ["NUM_SAMP"] + ["SAMP_CONSIS"]
    pandas = pd.DataFrame( index=list(clinic_data.index), columns=col )

    for patient in clinic_data.index:

        #access all samples of patient p
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

def add_specific_mutation_information( sample_data, patient_data, maf_data, output_file_location, HGVSp_Short="p.G12C", gene="KRAS" ):
    """Add specific mutation. TRUE/FALSE and counts per patient."""

    #check that all tumor samples in maf are also in sample data
    assert all( samp in list(sample_data.index) for samp in np.unique( maf_data["Tumor_Sample_Barcode"] )  )

    maf_data = maf_data[ maf_data["Hugo_Symbol"] == gene ]

    #find specific mutation
    hgvsp_mask = maf_data["HGVSp_Short"] == HGVSp_Short

    #tumor samples with mutation
    hgvsp_samples = list( maf_data["Tumor_Sample_Barcode"][ hgvsp_mask ] )

    #checking if code is same
    assert hgvsp_samples == list( maf_data[ hgvsp_mask ]["Tumor_Sample_Barcode"] )
    #check all samples are unique
    assert len( np.unique( hgvsp_samples) ) == len(hgvsp_samples)

    #create new column per SAMPLE to put number of mutation counts
    hgvsp_bool_col = "{}_PRESENT".format( HGVSp_Short )
    sample_data.loc[:, hgvsp_bool_col ] = [False] * sample_data.shape[0]
    sample_data.loc[ hgvsp_samples, hgvsp_bool_col ] = True
    assert sum( sample_data[hgvsp_bool_col] ) == len(hgvsp_samples)

    #map samples to patients. Only care about unique patients. If any sample of a patient is pos then pos
    patients = np.unique( sample_data["PATIENT_ID"][ sample_data[hgvsp_bool_col] ], return_counts=True )
    patient_data_mask_mut = []
    patient_data_mask_mut = [ True if p in patients[0] else False for p in patient_data.index ]
    assert len(patient_data_mask_mut) == patient_data.shape[0]
    patient_data[ hgvsp_bool_col ] = patient_data_mask_mut

    # inconsis = []
    # for p in patients[0][ patients[1] > 1 ]:
    #     patient_samples = sample_data[ sample_data["PATIENT_ID"] == p ]
    #     print(p)
    #     print(patient_samples)
    #     print(sum( patient_samples[hgvsp_bool_col] ))
    #     if sum( patient_samples[hgvsp_bool_col] ) == 0  or \
    #                                 sum(patient_samples[hgvsp_bool_col] ) == patient_samples.shape[0]:
    #         continue
    #     inconsis.append(p)
    #
    # consis_col_name = "{}_CONSIS".format( HGVSp_Short )
    # patient_data[consis_col_name] = [True] * patient_data.shape[0]
    # patient_data.loc[inconsis, consis_col_name] = False

    samp_output = "{}/genie_v8_kras_mutation_hgvs_sample.csv".format(output_file_location)
    sample_data.to_csv( samp_output )

    pat_output = "{}/genie_v8_kras_mutation_hgvs_patient.csv".format(output_file_location)
    patient_data.to_csv( pat_output )


def create_patientsXgene_count_pandas( mutation_data, sample_data, output_file ):
    """Create a patient by gene name mutation matrix. First column is KRASG12C. Mutations are counted"""

    assert all( samp in list(sample_data.index) for samp in np.unique( mutation_data["Tumor_Sample_Barcode"] )  )

    patient_list = [ sample_data.loc[sample]["PATIENT_ID"] for sample in mutation_data["Tumor_Sample_Barcode"] ]
    assert len(patient_list) == mutation_data.shape[0]

    mutation_data["PATIENT_ID"] = patient_list
    gene_list = list( np.unique( mutation_data["Hugo_Symbol"] ) )
    gene_list.append("KRAS.G12C")

    final_df = pd.DataFrame(0, index=np.unique(patient_list) , columns=gene_list )

    for patient in final_df.index:
        genes_mutated = list( mutation_data[ mutation_data[ "PATIENT_ID" ] == patient][ "Hugo_Symbol" ]  )
        if "KRAS" in genes_mutated:
            df_patient = mutation_data[ mutation_data["PATIENT_ID"] == patient ]
            df_kras = df_patient[ df_patient["Hugo_Symbol"] == "KRAS" ]
            final_df.loc[patient, "KRAS.G12C"] = sum( df_kras["HGVSp_Short"] == "p.G12C" )
        #print(len(genes_mutated))
        for gene in genes_mutated:
            final_df.loc[patient, gene] += 1

    csv_file = "{}/genie_v8_patients_X_genes_binary.csv".format(output_file)
    final_df.to_csv(csv_file)

if __name__ == '__main__':
    sd = pd.read_csv("GENIE_data/GENIE_v8/genie_v8_data_clinical_sample_colorectal.csv", index_col=0)
    #cd = pd.read_csv("GENIE_data/GENIE_v8/genie_v8_data_clinical_patient_colorectal.csv", index_col=0)
    maf = pd.read_csv("GENIE_data/GENIE_v8/genie_v8_data_mutations_extended_colorectal.csv", index_col=0, low_memory=False)
    create_patientsXgene_count_pandas(mutation_data = maf, sample_data = sd, output_file = "GENIE_data/GENIE_v8")
    # #create_anndata_mut_counts( sample_data=sd, clinic_data=cd ,maf_data=maf, \
    # #                           h5ad_output_file= "GENIE_data/GENIE_v8/genie_v8_anndata_mutation_count.h5ad")
    # #create_kras_mut_counts( sample_data=sd, clinic_data=cd ,maf_data=maf, \
    # #                             output_file= "GENIE_data/GENIE_v8/genie_v8_kras_mutation_count_sample.csv")
    # kras_count = pd.read_csv("GENIE_data/GENIE_v8/genie_v8_kras_mutation_count_sample.csv", index_col=0)
    # create_kras_mut_counts_per_patient( kras_mut_counts_sample = kras_count, clinic_data=cd, \
    #                             output_file = "GENIE_data/GENIE_v8/genie_v8_kras_mutation_count_patient.csv" )

    # sd = pd.read_csv("GENIE_data/GENIE_v8/genie_v8_kras_mutation_count_sample.csv", index_col=0)
    # cd = pd.read_csv("GENIE_data/GENIE_v8/genie_v8_kras_mutation_count_patient.csv", index_col=0)
    # add_specific_mutation_information( sd, cd, maf, "GENIE_data/GENIE_v8", HGVSp_Short="p.G12C" )
