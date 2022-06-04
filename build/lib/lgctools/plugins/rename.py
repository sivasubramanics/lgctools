import copy


def rename_data(smdata, msdata, samplemap_dict):
    out_msdata = copy.deepcopy(msdata)
    out_smdata = copy.deepcopy(smdata)

    for marker in msdata:
        for sample in msdata[marker].data:
            if sample in samplemap_dict:
                out_msdata[marker].data[samplemap_dict[sample]
                                        ] = out_msdata[marker].data.pop(sample)

    for sample in smdata:
        if sample in samplemap_dict:
            out_smdata[samplemap_dict[sample]] = out_smdata.pop(sample)

    return out_smdata, out_msdata
