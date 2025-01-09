#!/usr/bin/env python3

# ncbi-data.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import json
import csv
import re
import argparse

#----- ARGS -----#
parser = argparse.ArgumentParser()
parser.add_argument("--report", dest="report_file",  default="ncbi_dataset/data/data_report.jsonl", help="Path to 'data_report.jsonl' file from 'datasets download virus genome taxon'")
parser.add_argument("--subtype", dest="subtype_file",  default="ncbi-subtype.json", help="Path to JSON file from esearch'")
parser.add_argument("--taxids", dest="taxids_file",  default="ncbi-taxids.json", help="Path to JSON from 'datasets summary taxonomy taxon'")
parser.add_argument("--segsyns", dest="seg_syns", default=argparse.SUPPRESS,  help="Segment synonyms.")
args = parser.parse_args()

#----- FUNCTIONS -----#
# function for getting nested values from a dictionary
def get_nested_value(dictionary, keys, default='none'): 
    for key in keys: 
        dictionary = dictionary.get(key, default) 
        if dictionary == default: 
            break 
    return dictionary
# function for splitting segement synonyms into a dictionary
def splitSegSyns(segsyns):
    res = {}
    for seg in segsyns.split(';'):
        syns = seg.split('|')
        res[syns[0]] = list(set(syns + [ i.upper() for i in syns ] + [ i.lower() for i in syns ]))
    return res

# function for writing a dictionary to a CSV file
def write_dicts_to_csv(dicts, filename):
    # Get all unique keys
    keys = set()
    for d in dicts:
        keys.update(d.keys())
    keys = list(keys)

    # Write to CSV
    with open(filename, 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=keys)
        writer.writeheader()
        writer.writerows(dicts)

#----- MAIN -----#
# load NCBI Datasets taxid data
with open(args.taxids_file, 'r') as file:
    dataTaxa = json.load(file)
    speciesOptions = []
    for i in get_nested_value(dataTaxa, ['reports']):
        speciesOptions.append(get_nested_value(i, ['taxonomy','classification','species']))

# load NCBI Datasets report data
with open(args.report_file, 'r') as file:
    resReport = {}
    for line in file:
        dataReport = json.loads(line.strip())
        if 'accession' in dataReport:
            rowReport = {}
            rowReport['geographicRegion']   = get_nested_value(dataReport,['location','geographicRegion'])
            rowReport['organismName_host']  = get_nested_value(dataReport,['host','organismName'])
            rowReport['organismName_virus'] = get_nested_value(dataReport,['virus','organismName'])
            collectionDate                  = re.findall(r'\b\d{4}\b', get_nested_value(dataReport,['isolate','collectionDate']))
            rowReport['collectionDate']     = collectionDate[0] if collectionDate else 'none'
            taxIds                          = get_nested_value(dataReport,['virus','lineage'])
            if taxIds != 'none':
                for id in taxIds:
                    for option in speciesOptions:
                        if option['id'] == id['taxId']:
                            rowReport['species'] = option['name']
            else:
                rowReport['species'] = 'none'
            resReport[ dataReport['accession'] ] = rowReport

# load NCBI Esearch report data
## split segments, if supplied
segSyns = splitSegSyns(args.seg_syns) if 'seg_syns' in args else {}
with open(args.subtype_file, 'r') as file:
    resSubtype = {}
    segOpts    = []
    for line in file:
        dataSubtype = json.loads(line.strip())
        for i in dataSubtype.items():
            if 'result' in i:
                meta       = i[1]
                for uid in meta['uids']:
                    rowSubtype    = {}
                    subtypeKeys   = get_nested_value(meta[uid], ['subtype'])
                    subtypeValues = get_nested_value(meta[uid], ['subname'])
                    if subtypeKeys != 'none':
                        subtypeKeys   = subtypeKeys.split('|')
                        subtypeValues = subtypeValues.split('|')
                        counter       = 0
                        for key in subtypeKeys:
                            rowSubtype[key] = subtypeValues[counter]
                            counter         += 1
                    # update segment names based on supplied segment synonyms (if the are supplied)
                    if 'segment' in rowSubtype.keys():
                        segOpts.append(rowSubtype['segment'])
                        if len(segSyns) > 0:
                            for seg in segSyns.keys():
                                rowSubtype['segment'] = seg if rowSubtype['segment'] in segSyns[seg] else rowSubtype['segment']
                            rowSubtype['segment'] = rowSubtype['segment'] if rowSubtype['segment'] in segSyns.keys() else 'none'
                    else:
                        rowSubtype['segment'] = 'none'
                    resSubtype[ get_nested_value(meta[uid], ['accessionversion']) ] = rowSubtype
    print(f'Segment Options: {", ".join(set(segOpts))}')

# Combine datastreams
## all accessions in the NCBI Datasets report
resComplete = []
segComplete = []
targetKeys  = [ 'accession', 'geographicRegion', 'collectionDate', 'organismName_host', 'organismName_virus', 'species', 'subtype', 'genotype', 'serotype', 'lineage', 'clade', 'segment' ]
for acc in resReport.keys():
    if acc != 'none' and acc in resSubtype.keys():
        dictAll = resReport[acc] | resSubtype[acc] | {'accession': acc}
    else:
        dictAll = resReport[acc] | {'accession': acc}
    if 'segment' in dictAll.keys():
        segComplete.append(dictAll['segment'])
    else:
        dictAll.update({ 'segment': 'none' })
    resComplete.append({key: dictAll[key] for key in targetKeys if key in dictAll})
## determine if the organism is segmented & filter missing segment info
counter = 0
for row in resComplete:
    if len(set(segComplete)) > 0 and resComplete[counter]['segment'] != 'none':
        resComplete[counter] = resComplete[counter] | { 'segmented': 'true' }
    else:
        resComplete[counter]['segment'] = 'wg'
        resComplete[counter] = resComplete[counter] | { 'segmented': 'false' }
    counter              += 1
## accessions only in the NCBI Datasets report
resReportOnly = []
for acc in list(resReport.keys()):
    if acc not in list(resSubtype.keys()):
        resReportOnly.append( resReport[acc] | { 'accession': acc } )
## accessions only in the NCBI Esearch data
resSubtypeOnly = []
for acc in list(resSubtype.keys()):
    if acc not in list(resReport.keys()):
        resSubtypeOnly.append( resSubtype[acc] | { 'accession': acc } )

# Save to CSV file
write_dicts_to_csv(resComplete, 'ncbi-meta.complete.csv')
write_dicts_to_csv(resReportOnly, 'ncbi-meta.report-only.csv')
write_dicts_to_csv(resSubtypeOnly, 'ncbi-meta.subtype-only.csv')