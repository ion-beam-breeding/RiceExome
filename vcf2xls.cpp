/**
 * Copyright (c) 2017-2018 Hiroyuki Ichida. All rights reserved.
 *
 * @file vcf2xls.cpp
 * @author Hiroyuki Ichida <histfd@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 2 (GPL-2.0)
 * as published by the Free Software Foundation, Inc.
 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#include"histd.h"
#include"headerline.h"
#include<getopt.h>
#include<iostream>
#include<cstdlib>
#include<fstream>
#include<string>
#include<vector>
#include<map>
#include<set>
#include<climits>
#include<algorithm>

#define VCF2XLS_DEFAULT_COLUMNS     "GT,DP,AD"
#define VCF2XLS_SZ_FLANKS_SHOW      50
//------------------------------------------------------------------------------
typedef std::map<std::string,std::string> AnnotationDB;
typedef std::map<std::string,std::string> KeyValueDB;
typedef std::map<std::string,int> SampleOrderMap;
typedef std::vector<KeyValueDB> KeyValueDBArray;
typedef std::vector<hi::StringArray> StringArrayArray;
typedef std::vector<int> IntArray;
typedef std::map<std::string,std::string> MutEffectDB;
typedef std::map<int,std::string> GenotypeBaseMap;
typedef std::set<std::string> StringSet;
//------------------------------------------------------------------------------
inline bool get_record_as_string(const KeyValueDB &inputDB, const char *keyname, \
        std::string &result){

    KeyValueDB::const_iterator dbHit = inputDB.find(keyname);
    if(inputDB.end() != dbHit)
        result = dbHit->second;
    else{
        result = ".";
        return false;
    }

    return true;
}
//------------------------------------------------------------------------------
inline bool get_record_as_int(const KeyValueDB &inputDB, const char *keyname, \
        int *result){

    KeyValueDB::const_iterator dbHit = inputDB.find(keyname);
    if(inputDB.end() != dbHit)
        *result = std::atoi(dbHit->second.c_str());
    else{
        *result = -1;
        return false;
    }

    return true;
}
//------------------------------------------------------------------------------
bool parse_annotation(const std::string &line, KeyValueDB &result){
    hi::StringArray values;
    if(! hi::split(values, line, ';'))
        return false;

    std::string key, value;
    ssize_t sep;
    for(hi::StringArray::iterator iter=values.begin(); iter!=values.end(); ++iter){
        sep = iter->find('=');
        if(std::string::npos == sep){
            result.insert(std::pair<std::string,std::string>(*iter, *iter));
            continue;
        }
        key = iter->substr(0, sep);
        value = iter->substr(sep+1);
        result.insert(std::pair<std::string,std::string>(key, value));
    }

    return true;
}
//------------------------------------------------------------------------------
bool load_annotations_from_gff(const char *gff_fn, AnnotationDB &annots){

    std::ifstream infile(gff_fn, std::ios::in);
    if(infile.fail()){
        std::cerr << ERROR_STRING << "the specified GFF file (" << gff_fn \
                << ") can't open for reading." << ENDL;
        return false;
    }

    // process file
    std::string line, gene_id, gene_func;
    hi::StringArray items;
    KeyValueDB data;
    std::stringstream funcstr;
    while(std::getline(infile, line)){
        if('#' == line[0])
            continue;

        items.clear();
        hi::split(items, line, '\t');
        if(9 > items.size() || "gene" != items.at(2))
            continue;

        // parse
        data.clear();
        funcstr.str("");
        parse_annotation(items[8], data);

        // "ID"
        if(! get_record_as_string(data, "ID", gene_id))
            continue;

        // "CGSNL Gene Symbol"
        if(true==get_record_as_string(data, "CGSNL Gene Symbol", gene_func) \
                && "." != gene_func && "_" != gene_func)
            funcstr << "[[" << gene_func << "]] ";

        // "Note"
        if(! get_record_as_string(data, "Note", gene_func))
            gene_func = "NO_DESCRIPTION_AVAILABLE";

        funcstr << gene_func;

        // register
        annots.insert(std::pair<std::string,std::string>(gene_id, funcstr.str()));
    }

    infile.close();
    return true;
}
// -----------------------------------------------------------------------------
bool parse_info_field(const std::string &line, KeyValueDB &result){

    hi::StringArray values, arr;
    if(! hi::split(values, line, ';'))
        return false;

    std::string key, value;
    ssize_t sep;
    for(hi::StringArray::iterator iter=values.begin(); iter!=values.end(); ++iter){
        sep = iter->find('=');
        if(std::string::npos == sep){
            result.insert(std::pair<std::string,std::string>(*iter, *iter));
            continue;
        }
        key = iter->substr(0, sep);
        value = iter->substr(sep+1);
        result.insert(std::pair<std::string,std::string>(key, value));
    }

    return true;
}
//------------------------------------------------------------------------------
inline void write_output_header(const char *filename, const hi::StringArray &names, \
        const hi::StringArray &columns, const char *cmdstr, std::ostream &ofs){

    std::stringstream inputstr;
    inputstr << "input_fn=" << filename;
    write_basic_header(__FILE__, __DATE__, __TIME__, cmdstr, inputstr.str().c_str(), std::cout);
    ofs << "#CHROM" \
        << "\tChrStart" << "\tChrEnd" << "\tReference" << "\tAlternatives" \
        << "\tQuality" << "\tFilter" << "\tInfo" << "\tType" << "\tEffect" << "\tGene" \
        << "\tAnnotation";

    // Column titles for user-defined datasets
    for(hi::StringArray::const_iterator title=columns.begin(); title!=columns.end(); ++title){
        for(hi::StringArray::const_iterator name=names.begin(); name!=names.end(); ++name)
            ofs << '\t' << *name << '_' << *title;
    }
    ofs << "\tLink" << ENDL;
}
//------------------------------------------------------------------------------
inline std::string remove_string(const std::string &input, \
        const char *word, const int szWord){

    std::string target = input;
    ssize_t hit = target.find(word);
    while(std::string::npos != hit){
        target.replace(hit, szWord, "");
        hit = target.find(word);
    }

    return target;
}
//------------------------------------------------------------------------------
bool process_ann_field(const std::string &ann_string, \
        std::string &mut_effect, hi::StringArray &mut_genes){

    hi::StringArray blocks;
    hi::split(blocks, ann_string, ',');
    if(1 > blocks.size())
        return false;

    hi::StringArray items, parts;
    std::string geneid;
    for(hi::StringArray::iterator target=blocks.begin(); target!=blocks.end(); ++target){
        items.clear();
        hi::split(items, *target, '|');

        // effect
		if(1 < items.size() && "" != items[1])
            mut_effect = items[1];
        else{
	    	mut_effect = ".";
            return true;
	    }

        // mutated genes
		if(3 < items.size() && "" != items[3]){
			geneid = items[3];
        	geneid = remove_string(geneid, "LOC_", 4);
        	geneid = remove_string(geneid, "Gene_", 5);
        }
        else{
        	geneid = ".";
            break;
        }

        parts.clear();
        hi::split(parts, geneid, '-');
        if(4 < parts[0].length() && 't' == parts[0][4] && "Os" == parts[0].substr(0, 2))
            parts[0][4] = 'g';
        mut_genes.push_back(parts[0]);

        if(1 < parts.size() && 4 < parts[parts.size()-1].length() && \
                't' == parts[parts.size()-1][4] & "Os" == parts[parts.size()-1].substr(0, 2)){
            parts[parts.size()-1][4] = 'g';
            mut_genes.push_back(parts[parts.size()-1]);
        }

        // Only use the first hit -- use for() loop for future extension
        if(target==blocks.begin())
            break;
    }

    return true;
}
//------------------------------------------------------------------------------
bool process_ann_field(const std::string &ann_string, \
        std::string &effect_string, MutEffectDB &effectDB){

    hi::StringArray blocks;
    hi::split(blocks, ann_string, ',');
    if(1 > blocks.size())
        return false;

    hi::StringArray items, parts;
    std::string mut_effect;
    StringSet effects;
    for(hi::StringArray::iterator target=blocks.begin(); target!=blocks.end(); ++target){
        items.clear();
        hi::split(items, *target, '|');

        // effect
		if(1 < items.size() && "" != items[1]){
	        mut_effect = items[1];
            effects.insert(mut_effect);
        }
        else{
            effect_string = ".";
            return true;
        }

        // register only the first (most serious) effect for each variant
        // expect MutEffectDB is a std::map type
        effectDB.insert(std::pair<std::string,std::string>(items[0], mut_effect));
    }

    std::stringstream sstr;
    for(StringSet::const_iterator iter=effects.begin(); iter!=effects.end(); ++iter){
        if(effects.begin() != iter)
            sstr << ", " << *iter;
        else
            sstr << *iter;
    }
    effect_string = sstr.str();
    return true;
}
//------------------------------------------------------------------------------
bool create_genotype_base_map(const std::string &refString, const std::string &altString, \
        GenotypeBaseMap &genotypeBaseMap){

    // REF is always genotype '0'
    genotypeBaseMap.insert(std::pair<int,std::string>(0, refString));

    // ALT genotypes
    hi::StringArray arr;
    hi::split(arr, altString, ',');
    for(int gtid=0; gtid<arr.size(); ++gtid)
        genotypeBaseMap.insert(std::pair<int,std::string>(gtid+1, arr.at(gtid)));

    return true;
}
//------------------------------------------------------------------------------
inline void correct_mutation_type(const std::string &ref, const std::string &alt, \
        std::string &mut_type){

    if(1==ref.length() && 1==alt.length())
        mut_type = "SNV";
    else if(std::string::npos != alt.find(",") || std::string::npos != ref.find(","))
        mut_type = "MULTI_ALLELIC";
    else if(ref.length() > alt.length())
        mut_type = "DEL";
    else if(ref.length() < alt.length())
        mut_type = "INS";
    else if(ref.length() == alt.length())
        mut_type = "DEL+INS";
}
//------------------------------------------------------------------------------
inline void correct_mutation_length(const std::string &ref, const std::string &alt, \
        int *mut_length){

    if(std::string::npos != alt.find(",") || std::string::npos != ref.find(",")){
        hi::StringArray elements;
        hi::split(elements, alt, ',');
        hi::split(elements, ref, ',');
        int maxLength = 0;
        for(hi::StringArray::iterator iter=elements.begin(); iter!=elements.end(); ++iter)
            maxLength = std::max(maxLength, std::atoi(iter->c_str()));
        *mut_length = maxLength;
    }
    else
        *mut_length = alt.length() - ref.length();
}
//------------------------------------------------------------------------------
inline void get_annotation(const AnnotationDB &annotDB, const std::string &geneid, \
        std::string &result){

    AnnotationDB::const_iterator hit = annotDB.find(geneid);
    if(annotDB.end() != hit)
        result = hit->second;
    else
        result = ".";
}
//------------------------------------------------------------------------------
inline void generate_annotation_string(const AnnotationDB &annotDB, \
        const hi::StringArray &mut_genes, std::string &gidstr, std::string &funcstr){

    if(0 >= mut_genes.size()){
        gidstr  = ".";
        funcstr = ".";
        return;
    }

    gidstr = mut_genes.at(0);
    get_annotation(annotDB, gidstr, funcstr);

    if(2 <= mut_genes.size()){
        const std::string mutgene2 = mut_genes.at(mut_genes.size()-1);
        std::stringstream sstr;
        sstr << gidstr << " // " << mutgene2;
        gidstr = sstr.str();

        std::string func2;
        get_annotation(annotDB, mutgene2, func2);
        sstr.str("");
        sstr << funcstr << " // " << func2;
        funcstr = sstr.str();
    }
}
//------------------------------------------------------------------------------
inline void fetch_user_defined_datasets(const KeyValueDBArray &allSampleData, \
        const hi::StringArray &columnsToWrite, StringArrayArray &result){

    hi::StringArray dataset;
    KeyValueDB::const_iterator iter;
    for(KeyValueDBArray::const_iterator DB=allSampleData.begin(); DB!=allSampleData.end(); ++DB){
        dataset.clear();
        for(hi::StringArray::const_iterator key=columnsToWrite.begin(); key!=columnsToWrite.end(); ++key){
            iter = DB->find(*key);
            if(DB->end() != iter)
                dataset.push_back(iter->second);
            else
                dataset.push_back(".");
        }
        result.push_back(dataset);
    }
}
//------------------------------------------------------------------------------
inline void modify_GT(KeyValueDB &DB){

    // Switch
    bool isUseAD = true;

    // Retrive records
    KeyValueDB::iterator iterGT = DB.find("GT");
    if(DB.end() == iterGT)
        return;
    KeyValueDB::iterator iterAD = DB.find("AD");
    if(DB.end() == iterAD)
        isUseAD = false;

    // Parse AD items
    hi::StringArray ad_items;
    if(isUseAD)
        hi::split(ad_items, iterAD->second, ',');

    // Replace
    size_t hit = iterGT->second.find("/");
    while(std::string::npos != hit){
        iterGT->second.replace(hit, 1, "|");
        hit = iterGT->second.find("/");
    }

    // Correct genotypes
    if(".|." == iterGT->second){
        iterGT->second = ".";
        if(isUseAD){
            ad_items.clear();
            for(int i=0; i<2; ++i)
                ad_items.push_back(".");
        }
    }
    else if("1|." == iterGT->second){
        iterGT->second = ".";
        if(isUseAD){
            const std::string adstr = ad_items.at(0);
            ad_items.clear();
            ad_items.push_back(".");
            ad_items.push_back(adstr);
        }
    }

    // Genotype should be homozygous if the other AD is zero
    if(isUseAD && 2<=ad_items.size()){
        if("0"==ad_items.at(0) && "0"!=ad_items.at(1))
            iterGT->second = "1|1";
        else if("0"==ad_items.at(1) && "0"!=ad_items.at(0))
            iterGT->second = "0|0";
    }

    return;
}
//------------------------------------------------------------------------------
inline void correct_pindel_AD_format(KeyValueDB &DB, std::string &record){
    // Assume the record has single number in AD
    std::stringstream sstr;
    KeyValueDB::iterator gt = DB.find("GT");
    if(DB.end() == gt)
        return;

    if("0|0" == gt->second || "0/0" == gt->second){
        sstr << record << "|0";
        record = sstr.str();
    }
    else if("1|1" == gt->second || "1/1" == gt->second){
        sstr << "0|" << record;
        record = sstr.str();
    }
    else if("1|." == gt->second){
        gt->second = ".";
        sstr << ".|" << record;
        record = sstr.str();
    }

    return;
}
//------------------------------------------------------------------------------
inline void modify_AD(KeyValueDB &DB){

    KeyValueDB::iterator ad = DB.find("AD");
    if(DB.end() == ad)
        return;

    // genotype should always be ".(undetermined)" when DP == "."
    KeyValueDB::iterator dp = DB.find("DP");
    if(DB.end() != dp  &&  "." == dp->second){
        ad->second = ".";
        return;
    }

    ssize_t hit = ad->second.find(",");
    if(std::string::npos == hit){
        correct_pindel_AD_format(DB, ad->second);
        return;
    }

    while(std::string::npos != hit){
        ad->second.replace(hit, 1, "|");
        hit = ad->second.find(",");
    }
    return;
}
//------------------------------------------------------------------------------
inline void interpolate_DP(KeyValueDB &DB){

    // do nothing if DP is already in the record
    KeyValueDB::iterator iter = DB.find("DP");
    if(DB.end() != iter)
        return;

    // interpolate DP from AD (if avairable)
    iter = DB.find("AD");
    if(DB.end() == iter)
        return;

    hi::StringArray items;
    hi::split(items, iter->second, '|');

    int totalDepth = 0;
    for(hi::StringArray::iterator item=items.begin(); item!=items.end(); ++item)
        totalDepth += std::atoi(item->c_str());

    std::stringstream depthstr;
    depthstr << totalDepth;
    DB.insert(std::pair<std::string,std::string>("DP", depthstr.str()));
    return;
}
//------------------------------------------------------------------------------
inline void modify_data(KeyValueDBArray &allSampleData){
    for(KeyValueDBArray::iterator DB=allSampleData.begin(); DB!=allSampleData.end(); ++DB){
        // DO NOT change the order!!
        modify_GT(*DB);
        modify_AD(*DB);
        interpolate_DP(*DB);
    }
}
//------------------------------------------------------------------------------
bool parse_sample_fields(const hi::StringArray &items, const IntArray &sampleOrder, \
        KeyValueDBArray &result){

    // parse FORMAT field
    hi::StringArray format;
    hi::split(format, items.at(8), ':');
    if(0 >= format.size())
        return false;

    hi::StringArray elements;
    KeyValueDB DB;
    int minItems;
    for(IntArray::const_iterator column=sampleOrder.begin(); column!=sampleOrder.end(); ++column){
        elements.clear();
        hi::split(elements, items.at(*column), ':');
/*
        if(format.size() != elements.size()){
            std::cerr << WARNING_STRING << "invalid data structure:" \
                << " the number of elements differ between FORMAT and SAMPLE fields." \
                << ENDL;
        }
*/
        DB.clear();
        minItems = std::min(format.size(), elements.size());
        for(int k=0; k<minItems; ++k)
            DB.insert(std::pair<std::string,std::string>(format.at(k), elements.at(k)));
        result.push_back(DB);
    }
    return true;
}
//------------------------------------------------------------------------------
int get_min_fields(const StringArrayArray &input){

    int minFields = input.at(0).size();
    for(int i=1; i<input.size(); ++i)
        minFields = std::min(minFields, (int)input.at(i).size());
    return minFields;
}
//------------------------------------------------------------------------------
inline void trim_sequence(std::string &seq, const int *szFlanking){
    const size_t szSeq = seq.length();
    if((2 * (*szFlanking)) >= szSeq)
        return;

    std::stringstream sstr;
    sstr << seq.substr(0, *szFlanking) \
            << " ...(" << szSeq - (2 * (*szFlanking)) << " bp)... " \
            << seq.substr(szSeq-(*szFlanking));
    seq = sstr.str();
}
//------------------------------------------------------------------------------
inline void trim_annotation(std::string &Line, const char *Key){

    size_t headpos = Line.find(Key);
    if(std::string::npos == headpos)
        return;

    size_t endpos = Line.find(';', headpos);
    if(std::string::npos == endpos)
        Line.erase(headpos, std::string::npos);
    else
        Line.erase(headpos, endpos-headpos+1);

    if(';' == Line[Line.length()-1]){
        size_t pos = Line.find_last_not_of(";\n");
        if(std::string::npos != pos)
            Line.erase(pos+1, std::string::npos);
        else
            Line = ".";
    }
}
//------------------------------------------------------------------------------
void create_sample_order_array(const hi::StringArray &samples, \
        const SampleOrderMap &sampleOrderMap, IntArray &sampleOrder){

    for(hi::StringArray::const_iterator name=samples.begin(); name!=samples.end(); ++name){
        SampleOrderMap::const_iterator iter = sampleOrderMap.find(*name);
        if(iter!=sampleOrderMap.end())
            sampleOrder.push_back(iter->second);
    }
}
//------------------------------------------------------------------------------
bool process_vcf(const char *filename, const hi::StringArray &samples, \
        const hi::StringArray &columnsToWrite, const AnnotationDB &annotDB, \
        const int *szFlanking, const char *cmdstr){

    std::ifstream infile(filename, std::ios::in);
    if(infile.fail()){
        std::cerr << ERROR_STRING << "input file (" << filename << ") open failed." << ENDL;
        return false;
    }

    // process file
    IntArray sampleOrder;
    std::string line, chrName, ann, mut_type, mut_effect, mut_geneid, mut_genefunc;
    int mut_length, startPosition, endPosition;
    hi::StringArray items, mut_genes;
    KeyValueDB infoDB;
    KeyValueDBArray allSampleData;
    StringArrayArray dataToWrite;

    // write header
    write_output_header(filename, samples, columnsToWrite, cmdstr, std::cout);

    while(std::getline(infile, line)){
        if("##" == line.substr(0, 2))
            continue;
        else if("#CHROM" == line.substr(0, 6)){
            hi::StringArray elements;
            hi::split(elements, line, '\t');
            if(9 > elements.size()){
                std::cerr << "invalid header structure. line=" << line << ENDL;
                return false;
            }

            sampleOrder.clear();
            sampleOrder.resize(samples.size());
            for(size_t i=9; i<elements.size(); ++i){
                hi::StringArray::const_iterator hit = std::find( \
                            samples.begin(), samples.end(), elements.at(i));
                if(samples.end() != hit){
                    int pos = std::distance(samples.begin(), hit);
                    sampleOrder.at(pos) = i;
                }
            }
            continue;
        }
        items.clear();
        hi::split(items, line, '\t');
        if(7 > items.size()){
            std::cerr << WARNING_STRING << ENDL;
            continue;
        }

        // parse INFO field
        infoDB.clear();
        parse_info_field(items[7], infoDB);

        // SVTYPE & SVLEN
        get_record_as_string(infoDB, "SVTYPE", mut_type);
        if("." == mut_type)
            correct_mutation_type(items[3], items[4], mut_type);
        get_record_as_int(infoDB, "SVLEN", &mut_length);
        if(0 > mut_length)
            correct_mutation_length(items[3], items[4], &mut_length);

        // Location
        chrName = items[0];
        startPosition = std::atoi(items[1].c_str());
        endPosition = startPosition + std::abs(mut_length);

        // ANN
        mut_genes.clear();
        mut_effect = ".";
        if(get_record_as_string(infoDB, "ANN", ann))
            process_ann_field(ann, mut_effect, mut_genes);

        // REF/ALT sequences
        trim_sequence(items[3], szFlanking);
        trim_sequence(items[4], szFlanking);

        // Remove some features from INFO
        trim_annotation(items[7], "ANN");
        trim_annotation(items[7], "LOF");

        // Output
        // CHROM, StartPos, EndPos, REF, ALT, QUAL, FILTER, INFO, Type, Effect
        std::cout << chrName << '\t' << startPosition << '\t' << endPosition \
                << '\t' << items[3] << '\t' << items[4] << '\t' << items[5] \
                << '\t' << items[6] << '\t' << items[7] << '\t' << mut_type \
                << '\t' << mut_effect;

        // Gene & Annotation
        generate_annotation_string(annotDB, mut_genes, mut_geneid, mut_genefunc);
        std::cout << '\t' << mut_geneid << '\t' << mut_genefunc;

        // Create FORMAT->DATA map from each sample
        allSampleData.clear();
        parse_sample_fields(items, sampleOrder, allSampleData);
        modify_data(allSampleData);

        // Extract user-defined datasets
        dataToWrite.clear();
        fetch_user_defined_datasets(allSampleData, columnsToWrite, dataToWrite);
        for(int i=0; i<get_min_fields(dataToWrite); ++i){
            for(StringArrayArray::const_iterator fields=dataToWrite.begin(); fields!=dataToWrite.end(); ++fields)
                std::cout << '\t' << fields->at(i);
        }

        // Link
        std::cout << '\t' << "=HYPERLINK(\"http://localhost:60151/goto?locus=" \
            << chrName << ':' << startPosition << '-' << endPosition \
            << "\", \"link\")" << ENDL;
    }

    infile.close();
    return true;
}
//------------------------------------------------------------------------------
inline void print_usage(const char *cmd){
    std::cerr << "usage:" << ENDL;
    std::cerr << cmd \
            << " -i (vcf_fn) -a (annotation_gff) " \
            << " -c (column_to_export) -s (sz_flanking_to_show)" \
            << " [sample names]" << ENDL;
    std::cerr << ENDL;
}
// -----------------------------------------------------------------------------
bool determine_sample_names(const char *vcf_fn, hi::StringArray &names){

    std::ifstream infile(vcf_fn, std::ios::in);
    if(infile.fail()){
        std::cerr << ERROR_STRING << "input file (" << vcf_fn << ") open failed." << ENDL;
        return false;
    }

    std::string line;
    names.clear();
    while(std::getline(infile, line)){
        if("#CHROM" == line.substr(0, 6)){
            hi::StringArray items;
            hi::split(items, line, '\t');
            if(9 > items.size()){
                std::cerr << ERROR_STRING << "invalid header structure. line=" << line << ENDL;
                return false;
            }
            names.insert(names.end(), items.begin()+9, items.end());
            break;
        }
    }

    if(names.empty())
        return false;
    return true;
}
// -----------------------------------------------------------------------------
int main(int argc, char *argv[]){

    // parse arguments
    int nArg=2;
    if(argc<nArg){
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    std::string gff_fn="", vcf_fn="", columnStr=VCF2XLS_DEFAULT_COLUMNS;

    // parse arguments
    char option;
    int szFlanking=VCF2XLS_SZ_FLANKS_SHOW;
    while ((option = getopt(argc, argv, "i:a:c:s:")) != -1){
        switch (option){
            case 'i':
                vcf_fn = optarg;
                break;
            case 'a':
                gff_fn = optarg;
                break;
            case 'c':
                columnStr = optarg;
                break;
            case 's':
                szFlanking = std::atoi(optarg);
                break;
            default:
                print_usage(argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    // vcf_fn must be specified
    if("" == vcf_fn){
        std::cerr << ERROR_STRING << "input VCF (-i) must be specified." << ENDL;
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    // szFlanking must be 1 or larger value
    if(0 >= szFlanking){
        szFlanking = VCF2XLS_SZ_FLANKS_SHOW;
        std::cerr << WARNING_STRING << "invalid szFlanking size specified." \
                << "Using the default value (" << VCF2XLS_SZ_FLANKS_SHOW << ")." << ENDL;
    }

    // column order
    hi::StringArray names;
    SampleOrderMap sampleOrderMap;
    for(int i=optind; i<argc; ++i){
        names.push_back(argv[i]);
        sampleOrderMap.insert(std::pair<std::string,int>(argv[i], i-optind));
    }
    // automatically determine sample names if not specified
    if(names.empty()){
        std::cerr << INFO_STRING << "Sample names not specified." \
                << " Retriving from the header line...";
        if(determine_sample_names(vcf_fn.c_str(), names)){
            std::sort(names.begin(), names.end(), std::less<std::string>());
            std::cerr << " Found!" << ENDL;
            std::cerr << INFO_STRING << "Followings are subject to process:";
            for(hi::StringArray::iterator iter=names.begin(); iter!=names.end(); ++iter)
                std::cerr << ' ' << *iter;
            std::cerr << ENDL;
        }
        else
            std::cerr << " NOT FOUND.";
    }

    // user-defined datasets to export
    hi::StringArray columns;
    if("" != columnStr)
        hi::split(columns, columnStr, ',');

    // load annotatinos from GFF
    AnnotationDB annots;
    if("" != gff_fn)
        load_annotations_from_gff(gff_fn.c_str(), annots);

    // process file
    std::string cmdstr = generate_cmd_string(argc, argv);
    process_vcf(vcf_fn.c_str(), names, columns, annots, &szFlanking, cmdstr.c_str());

    exit(EXIT_SUCCESS);
}
// -----------------------------------------------------------------------------
