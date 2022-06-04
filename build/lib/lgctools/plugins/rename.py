
def rename_data(smdata, msdata, samplemap_dict):
    out_msdata = msdata.copy()
    out_smdata = smdata.copy()

    for marker in msdata:
        for sample in msdata[marker].data:
            print(f"{sample}\t{marker}")
            if sample in samplemap_dict:
                out_msdata[marker].data[samplemap_dict[sample]
                                        ] = out_msdata[marker].data.pop(sample)

    for sample in smdata:
        if sample in samplemap_dict:
            out_smdata[samplemap_dict[sample]] = out_smdata.pop(sample)

    return out_smdata, out_msdata
