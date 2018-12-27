/**
 * Copyright (c) 2015-2018 Hiroyuki Ichida. All rights reserved.
 *
 * @file genotype_filter.h
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

#ifndef EXOME_GENOTYPE_FILTER_H
#define EXOME_GENOTYPE_FILTER_H

#include "histd.h"
#include "headerline.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdio>
#include <map>
#include <set>

#define KEY_00      "0|0"
#define KEY_11      "1|1"
#define KEY_NULL    "."     // vcf2xlsを使用して作製したファイルの場合
#define KEY_EMPTY   ""      // オリジナルの理研ジェネシス書式の場合
#define KEY_GT      "_GT"
#define KEY_DP      "_DP"
#define KEY_AD      "_AD"
#define KEY_CHROM       "CHROM"
#define KEY_PROGRAM     "Program"
#define DEFAULT_PROGRAM_NAME    "Program"
#define NUM_HETERO_ALLOWED  0
#define NUM_NULL_ALLOWED    0
#define MIN_SUPPORT_READS   10
#define MAX_ALT_RATE        0.05
#define MAX_ALT_LINE_FRAC   0.5
#define MAX_HETERO_BIAS     3.0     // float
//------------------------------------------------------------------------------
typedef std::vector<int> IntArray;
typedef std::set<std::string> StringSet;
typedef std::map<int, std::string> ColumnToName;
//------------------------------------------------------------------------------
/**
 * 遺伝子型がホモ接合か否かを検定する
 * @param genotype 遺伝子型（3文字，区切り文字を挟んだ両端を比較する）
 * @return 遺伝子型がホモ接合か否か
 */
bool is_homozygous(const std::string &genotype){

    if(genotype[0]==genotype[2])
        return true;

    return false;
}
//------------------------------------------------------------------------------
/**
 * 遺伝子型が同じ対立遺伝子のみからなるか否かを検定する
 * @param gtype1 遺伝子型1（3文字，区切り文字を挟んだ両端を比較する）
 * @param gtype2 遺伝子型2（3文字，区切り文字を挟んだ両端を比較する）
 * @return 与えられた2つの遺伝子型が同じ対立遺伝子のみからなるか否か
 */
bool is_by_same_alleles(const std::string &gtype1, const std::string &gtype2){

    if('0'!=gtype1[0] && gtype1[0]!=gtype2[0] && gtype1[0]!=gtype2[2])
        return false;
    if('0'!=gtype1[2] && gtype1[2]!=gtype2[0] && gtype1[2]!=gtype2[2])
        return false;

    return true;
}
//------------------------------------------------------------------------------
/**
 * 変異が系統特異的か否かを検定する
 * @param  arr                 行をタブ区切りでパーズした配列
 * @param  gtColumns           判定に用いるGenotypeのカラムを列挙した配列
 * @param  max_hetero          許容するヘテロ接合の系統数
 * @param  null_genotype       Genotypeが不明と見なす文字列
 * @param  specificColumnIndex 特異的な系統のカラム位置
 * @return                     変異が系統特異的であるか否か
 */
bool is_line_specific(const hi::StringArray &arr, \
        const IntArray &gtColumns,  const int *max_hetero, \
        const int *max_null_allowed, const char *null_genotype, \
        int *specificColumnIndex){

    int mutant_pos_index=-1, num_hetero=0, num_null=0;
    bool isTargetHeterozygous=false, isHomo, isBySameAlleles;

    for(IntArray::const_iterator \
            pos=gtColumns.begin(); pos!=gtColumns.end(); ++pos){

        if(KEY_NULL==arr[*pos] || KEY_EMPTY==arr[*pos]){
            ++num_null;
            continue;
        }
        else if(null_genotype==arr[*pos])
            continue;

        if(*max_null_allowed < num_null){
            *specificColumnIndex = -1;
            return false;
        }

        isHomo = is_homozygous(arr[*pos]);
        if(0 > mutant_pos_index){
            // this means a specific mutant hasn't been found yet
            mutant_pos_index = std::distance(gtColumns.begin(), pos);
            if(! isHomo)
                isTargetHeterozygous = true;
            continue;
        }

        isBySameAlleles = is_by_same_alleles( \
                arr[gtColumns.at(mutant_pos_index)], arr[*pos]);

        if(isTargetHeterozygous && isHomo){
            // ヘテロ接合の個体における変異として登録されていて，
            // 後からホモ接合の個体が見つかった場合
            mutant_pos_index = std::distance(gtColumns.begin(), pos);
            isTargetHeterozygous = false;

            if(isBySameAlleles){
                ++num_hetero;
                if(num_hetero > *max_hetero){
                    *specificColumnIndex = -1;
                    return false;
                }
            }
            continue;
        }

        if(isHomo && isBySameAlleles){
            // 同じアリルの組み合わせからなるホモ接合体が見つかった場合
            // 無条件にreject
            *specificColumnIndex = -1;
            return false;
        }
        else if(isBySameAlleles){
            // 同じアリルのヘテロ接合体の場合
            // これまでに出現したヘテロの数がmax_heteroを超える場合はreject
            ++num_hetero;
            if(num_hetero > *max_hetero){
                *specificColumnIndex = -1;
                return false;
            }
        }
    }

    if(*max_null_allowed < num_null){
        *specificColumnIndex = -1;
        return false;
    }

    *specificColumnIndex = mutant_pos_index;
    return true;
}
//------------------------------------------------------------------------------
/**
 * keywordを含むカラムの位置を列挙する
 * @param  header_elements ヘッダー行をタブ区切りで分割した配列
 * @param  keyword         検索する単語
 * @param  ignoreNames     解析対象としない系統名を列挙した配列
 * @param  columnPos       keywordを末尾に含むカラムの位置
 * @return                 keywordを末尾に含むカラムが見つかったか否か
 */
bool find_column_positions(hi::StringArray &header_elements, \
        const char *keyword, const StringSet &ignoreNames, IntArray &columnPos){

    const int szKeyword = std::strlen(keyword);
    for(size_t pos=0; pos<header_elements.size(); pos++){
        if(0 == std::strncmp( \
                header_elements[pos].substr( \
                    header_elements[pos].length()-szKeyword, szKeyword).c_str(), \
                    keyword, szKeyword)){

            std::string name = \
                    header_elements[pos].substr( \
                        0, header_elements[pos].length()-szKeyword);
            if(ignoreNames.end() == ignoreNames.find(name))
                columnPos.push_back(pos);
        }
    }

    if(0 == columnPos.size())
        return false;
    return true;
}
//------------------------------------------------------------------------------
bool find_column_position(hi::StringArray &header_elements, \
        const std::string &keyword, int *column_pos){

    *column_pos = -1;
    ssize_t hit;
    for(size_t pos=0; pos<header_elements.size(); pos++){
    	hit = header_elements.at(pos).find(keyword);
    	if(std::string::npos != hit){
    	    *column_pos = pos;
    	    return true;
    	}
    }
    return false;
}
//------------------------------------------------------------------------------
/**
 * 対立遺伝子の数を調べる
 * @param adepth AllelicDepthの文字列
 * @return 対立遺伝子の数
 */
int get_num_alleles(const std::string &adepth){
    const char *adcstr = adepth.c_str();
    int count=1;
    for(int i=0; i<std::strlen(adcstr); ++i){
        if('|' == adcstr[i])
            ++count;
    }
    return count;
}
//------------------------------------------------------------------------------
/**
 * 行全体での対立遺伝子の最大数を求める
 * @param  elements  タブ区切りで行を分解した配列
 * @param  adColumns 調べるカラム（AllelicDepth: "_AD"）の位置
 * @return レコード全体の対立遺伝子の
 */
int get_max_alleles(const hi::StringArray &elements, const IntArray &adColumns){
    int num_alleles = 0;
    for(IntArray::const_iterator \
            pos=adColumns.begin(); pos!=adColumns.end(); ++pos){

        num_alleles = std::max(num_alleles, get_num_alleles(elements.at(*pos)));
    }

    return num_alleles;
 }
 //------------------------------------------------------------------------------
/**
 * ホモ接合にも関わらず検出される対立遺伝子の割合でフィルタリングする
 * @param adepth AllelicDepthの文字列
 * @return MAX_ALT_RATEを超えるリードが存在するか否か
 */
bool is_pass_alt_rate_filter(const std::string &adepth){

    hi::StringArray elements;
    hi::split(elements, adepth, '|');
    if(2 < elements.size())
        return false;

    int adepthL = std::atoi(elements[0].c_str());
    int adepthR = std::atoi(elements[1].c_str());

    float adepth1 = (float)std::max(adepthL, adepthR);
    float adepth2 = (float)std::min(adepthL, adepthR);

    if(adepth2 / adepth1 > MAX_ALT_RATE)
        return false;

    return true;
}
//------------------------------------------------------------------------------
/**
 * ヘテロ接合でリードの比率が1:1から著しく外れたレコードをフィルタリングする
 * @param adepth AllelicDepthの文字列
 * @return MAX_HETERO_BIASを超えるリードが存在するか否か
 */
bool is_pass_hetero_bias_filter(const std::string &adepth){

    hi::StringArray elements;
    hi::split(elements, adepth, '|');
    if(2 < elements.size())
        return false;

    int adepthL = std::atoi(elements[0].c_str());
    int adepthR = std::atoi(elements[1].c_str());
    float ratio = (float)std::max(adepthL, adepthR) \
            / ((float)std::min(adepthL, adepthR)+0.001);

    if(MAX_HETERO_BIAS < ratio)
        return false;

    return true;
}
//------------------------------------------------------------------------------
/**
 * 系統特異性を判定する（対立遺伝子が2種類のレコードのみを認めるモード)
 * @param elements           タブ区切りで行をパーズした配列
 * @param line               （出力用）パーズ前の行データ
 * @param names              Index->系統名の変換テーブル
 * @param program_name       （出力用）プログラム名
 * @param gtColumns          判定対象のGenotypeカラムの位置
 * @param adColumns          判定対象のAllelicDepthカラムの位置
 * @param nHeteroAllowed     許容するヘテロ接合の系統数
 * @param useAltRateFilter   Alt rate filterの適用有無
 * @param maxAltLineFrac     [description]
 * @param useReadDepthFilter Read depth filterの適用有無
 * @param minSupportReads    [description]
 * @param isInvertSelection  系統特異的でないレコードのみを出力する
 */
void proc_by_two_alleles_mode(const hi::StringArray &elements, \
        ColumnToName &names, const std::string &program_name, \
        const IntArray &gtColumns, const IntArray &adColumns, \
        const int *nHeteroAllowed, const int *nNullAllowed, \
        const bool *useAltRateFilter, const float *maxAltLineFrac, \
        const bool *useReadDepthFilter, const int *minSupportReads, \
        const bool *isInvertSelection){

    // combine items in the parsed array
    std::stringstream linestr;
    for(hi::StringArray::const_iterator \
            iter=elements.begin(); iter!=elements.end(); ++iter){

        linestr << '\t' << *iter;
    }
    const std::string line = linestr.str();

    // check line specificity
    int target_allele_pos=1, specific_column_index;
    is_line_specific(elements, gtColumns, nHeteroAllowed, nNullAllowed, \
            KEY_00, &specific_column_index);
    if(0 > specific_column_index){
        // a line had 0|0 or 0|1 and all others had 1|1
        is_line_specific(elements, gtColumns, nHeteroAllowed, nNullAllowed, \
                KEY_11, &specific_column_index);
        target_allele_pos = 0;
        if(0 > specific_column_index || KEY_11!=elements[specific_column_index]){
            if(*isInvertSelection)
                std::cout << program_name << "\t." << line << ENDL;
            return;
        }
    }
    if(*isInvertSelection)
        return;

    // alt rate filter
    if(*useAltRateFilter){
        int alt_rate_exceed=0, hetero_bias_exceed=0;
        int num_homozygous=0;
        for(int i=0; i<gtColumns.size(); ++i){
            const bool isHomozygous = is_homozygous(elements[gtColumns.at(i)]);
            if(i==specific_column_index){
                if(false==isHomozygous \
                        && false==is_pass_hetero_bias_filter(elements[adColumns.at(i)])){
                    return;
                }
                continue;
            }
            else if(isHomozygous){
                ++num_homozygous;
                if(! is_pass_alt_rate_filter(elements[adColumns.at(i)]))
                    ++alt_rate_exceed;
            }
        }
        if(0 < ((float)alt_rate_exceed/(float)num_homozygous) - *maxAltLineFrac)
            return;
    }

    // filter by the read depth of a mutant allele
    if(*useReadDepthFilter){
        hi::StringArray adepth_elements;
        hi::split(adepth_elements, elements[adColumns.at(specific_column_index)], '|');
        if(target_allele_pos >= adepth_elements.size() \
                || std::atoi(adepth_elements.at(target_allele_pos).c_str()) < *minSupportReads)
            return;
    }

    // output
    ColumnToName::iterator iter = names.find(specific_column_index);
    if(names.end()!=iter){
        std::cout << program_name << '\t' << iter->second << line << ENDL;
    }
    else{
        std::cerr << ERROR_STRING << "specific_pos=" << specific_column_index \
                << ", record=" << program_name << '\t' << iter->second << line << ENDL;
    }

    return;
}
//------------------------------------------------------------------------------
/**
 * 対立遺伝子が3種類以上のレコードを認めるモード（他品種等のリファレンス配列にマッピングした結果に使用）
 * @method proc_by_multiple_alleles_mode
 * @param  file           std::ifstream &file  **ヘッダーは予め飛ばしておくこと**
 * @param  allele_start   alleleカラムの開始位置(左端=0)
 * @param  allele_end     alleleカラムの終了位置の次
 * @return                処理が成功したか否か
 */
void proc_by_multiple_alleles_mode(const hi::StringArray &arr, \
        ColumnToName &names, const std::string &program_name, \
        const IntArray &gtColumns, const IntArray &adColumns, \
        const int *nHeteroAllowed, const int *nNullAllowed, \
        const bool *isInvertSelection){

    int count, num_homo, num_hetero, specific_column_index, specific_gtype;
    int num_alleles = get_max_alleles(arr, adColumns);

    // mark homozygous columns in boolean
    // all others are heterozygous
    bool *isHomo = new bool [gtColumns.size()];
    for(int i=0; i<gtColumns.size(); ++i){
        // ignore a record with an empty genotype
        // this is usually a result of insufficient depth
        if(KEY_NULL==arr[i] || KEY_EMPTY==arr[i])
            return;
        isHomo[i] = is_homozygous(arr[gtColumns.at(i)]);
    }

    // combine items in the parsed array
    std::stringstream linestr;
    for(hi::StringArray::const_iterator \
            iter=arr.begin(); iter!=arr.end(); ++iter){

        linestr << '\t' << *iter;
    }
    const std::string line = linestr.str();

    // test all possible genotypes
    char gtchr;
    hi::StringArray adepth_elements;
    for(int gtype=1; gtype<num_alleles; ++gtype){
        gtchr = char(gtype + 0x30);
        count = num_homo = num_hetero = 0;
        specific_column_index = -1;

        for(int i=0; i<gtColumns.size(); ++i){
            if(gtchr==arr[gtColumns.at(i)][2] || gtchr==arr[gtColumns.at(i)][0]){
                if(0==num_homo){
                    specific_column_index = i;
                    specific_gtype = gtype;
                }

                if(isHomo[i])
                    ++num_homo;
                else
                    ++num_hetero;
            }
        }

        // final decision
        if(0<=specific_column_index &&  \
            ((1==num_homo && *nHeteroAllowed>=num_hetero) \
                    || (0==num_homo && 1==num_hetero))){
            // filtering by alt rate
            /** IMPLEMENT HERE **/

            // filter by read depth of the mutant allele
            adepth_elements.clear();
            hi::split(adepth_elements, arr[adColumns.at(specific_column_index)], '|');
            if(MIN_SUPPORT_READS > std::atoi(adepth_elements.at(specific_gtype).c_str()))
                continue;

            // output
            ColumnToName::iterator iter = names.find(specific_column_index);
            if(names.end()!=iter){
                std::cout << program_name << '\t' << iter->second \
                        <<'\t' << line << ENDL;
            }
            else{
                std::cerr << ERROR_STRING << "specific_pos=" \
                    << specific_column_index << ", record="  \
                    << program_name << '\t' << iter->second <<'\t' << line << ENDL;
            }
        }
        else if(*isInvertSelection && 0 > specific_column_index)
            std::cout << program_name << "\t.\t" << line << ENDL;
    }

    return;
}
// -----------------------------------------------------------------------------
bool create_column_to_name(const hi::StringArray &headerElements, \
        const IntArray &gtColumns, const char *keyword, ColumnToName &names){

    const int szKeyword = std::strlen(keyword);
    for(IntArray::const_iterator \
            pos=gtColumns.begin(); pos!=gtColumns.end(); ++pos){

        std::string name = \
                headerElements[*pos].substr(0, headerElements[*pos].length()-szKeyword);
        names.insert(std::pair<int,std::string>(std::distance(gtColumns.begin(), pos), name));
    }
    return true;
}
// -----------------------------------------------------------------------------

#endif
