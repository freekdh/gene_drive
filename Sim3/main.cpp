/* Freek de Haas, UBC, freekdh@gmail.com */

#include <iostream>
#include <vector>
#include "random.h"
#include "utils.h"
#include <assert.h>
#include <fstream>
#include "recomb.h"

const int n_genotypes = 256;
void print(std::vector<int> & , std::ofstream &);
void print_header(std::ofstream &);
std::vector<int> initialize_pop(const int &, const int &);
rnd::discrete_distribution make_dist(const std::vector<int> &);
const std::vector<double>  initialize_multiplicative_payload(const double &sd, const double &spl, const double &e);
const std::vector<double> calculate_frequencies(const std::vector<int> &pop,const int &total, const double &d, const double &r);

/* main */
int main(int argc, char *argv[]){
    const int n_wild = atoi(argv[1]);
    const int n_drive = atoi(argv[2]);
    const int max_pop_size = atoi(argv[3]);
    const int n_gen = atoi(argv[4]);
    const double avg_offspring =atof(argv[5]);
    const double drive_rate = atof(argv[6]);
    const double recomb_rate = atof(argv[7]);
    const double selection_drive = atof(argv[8]);
    const double selection_payload = atof(argv[9]);
    const double selection_toxin = atof(argv[10]);

    rnd::set_seed();
    std::ofstream output_file("./Data/sim_results.csv");
    output_file.fill(',');
    print_header(output_file);
    const std::vector<double> type_payload = initialize_multiplicative_payload(selection_drive,selection_payload,selection_toxin);

    std::vector<int> parents = initialize_pop(n_wild, n_drive);
    print(parents,output_file);
    for(int gen = 0; gen < n_gen; ++gen){
        
        std::vector<int> offspring_after_selection(n_genotypes,0);
        int total = 0;
        for(int i = 0; i<n_genotypes; ++i){
            offspring_after_selection[i]=rnd::poisson(parents[i]*avg_offspring*type_payload[i]);
            total+= offspring_after_selection[i];
        }
        rnd::discrete_distribution dist = make_dist(offspring_after_selection);
        while(total > max_pop_size){
            const int remove_individual = dist.sample();
            --offspring_after_selection[remove_individual];
            dist[remove_individual] = dist[remove_individual]-1;
            --total;
        }

        std::vector<int> offspring_after_mating(n_genotypes,0);
        std::vector<double> freq_after_mating = calculate_frequencies(offspring_after_selection,total,drive_rate,recomb_rate);
        rnd::discrete_distribution dist2(n_genotypes);
        for(int i = 0; i<n_genotypes;++i){
            dist2[i] = freq_after_mating[i];
        }
        for(int i = 0; i < total; ++i){
            ++offspring_after_mating[dist2.sample()];
        }
        parents = offspring_after_mating;

        print(parents,output_file);
    }
    return 0;
}

const std::vector<double> calculate_frequencies(const std::vector<int> &pop,const int &total, const double &d, const double &r){
    std::vector<double> F(n_genotypes,0.0);
    for(int i = 0; i < n_genotypes; ++i){
        F[i] = pop[i]/(double)total;
    }

    static const double OneMinusdd = 1.0-d*d;
    static const double OneMinusddd = 1.0-d*d*d;
    static const double dd = d*d;
    static const double ddd = d*d*d;
    
    const std::vector<double> F2{
        F[0],F[1],F[2],F[3],F[4],-((-1 + d)*F[5]),-((-1 + d)*F[6]),OneMinusdd*F[7],F[8],F[9],F[10],F[11],-((-1 + d)*F[12]),OneMinusdd*F[13],OneMinusdd*F[14],OneMinusddd*F[15],F[16],F[17],F[18],F[19],
        -((-1 + d)*F[20]),d*(F[5] + F[20]) + F[21],OneMinusdd*F[22],-((-1 + d)*(d*(F[7] + F[22]) + F[23])),F[24],F[25],F[26],F[27],OneMinusdd*F[28],-((-1 + d)*(d*(F[13] + F[28]) + F[29])),OneMinusddd*F[30],
        OneMinusdd*(d*(F[15] + F[30]) + F[31]),F[32],F[33],F[34],F[35],-((-1 + d)*F[36]),OneMinusdd*F[37],d*(F[6] + F[36]) + F[38],-((-1 + d)*(d*(F[7] + F[37]) + F[39])),F[40],F[41],F[42],F[43],
        OneMinusdd*F[44],OneMinusddd*F[45],-((-1 + d)*(d*(F[14] + F[44]) + F[46])),OneMinusdd*(d*(F[15] + F[45]) + F[47]),F[48],F[49],F[50],F[51],OneMinusdd*F[52],-((-1 + d)*(d*(F[37] + F[52]) + F[53])),
        -((-1 + d)*(d*(F[22] + F[52]) + F[54])),dd*(F[7] + F[22] + F[37] + F[52]) + d*(F[23] + F[39] + F[53] + F[54]) + F[55],F[56],F[57],F[58],F[59],OneMinusddd*F[60],OneMinusdd*(d*(F[45] + F[60]) + F[61]),
        OneMinusdd*(d*(F[30] + F[60]) + F[62]),-((-1 + d)*(dd*(F[15] + F[30] + F[45] + F[60]) + d*(F[31] + F[47] + F[61] + F[62]) + F[63])),F[64],-((-1 + d)*F[65]),-((-1 + d)*F[66]),OneMinusdd*F[67],F[68],
        -((-1 + d)*F[69]),-((-1 + d)*F[70]),OneMinusdd*F[71],-((-1 + d)*F[72]),OneMinusdd*F[73],OneMinusdd*F[74],OneMinusddd*F[75],d*(F[12] + F[72]) + F[76],-((-1 + d)*(d*(F[13] + F[73]) + F[77])),
        -((-1 + d)*(d*(F[14] + F[74]) + F[78])),OneMinusdd*(d*(F[15] + F[75]) + F[79]),-((-1 + d)*F[80]),d*(F[65] + F[80]) + F[81],OneMinusdd*F[82],-((-1 + d)*(d*(F[67] + F[82]) + F[83])),-((-1 + d)*F[84]),
        d*(F[69] + F[84]) + F[85],OneMinusdd*F[86],-((-1 + d)*(d*(F[71] + F[86]) + F[87])),OneMinusdd*F[88],-((-1 + d)*(d*(F[73] + F[88]) + F[89])),OneMinusddd*F[90],OneMinusdd*(d*(F[75] + F[90]) + F[91]),
        -((-1 + d)*(d*(F[28] + F[88]) + F[92])),dd*(F[13] + F[28] + F[73] + F[88]) + d*(F[29] + F[77] + F[89] + F[92]) + F[93],OneMinusdd*(d*(F[30] + F[90]) + F[94]),
        -((-1 + d)*(dd*(F[15] + F[30] + F[75] + F[90]) + d*(F[31] + F[79] + F[91] + F[94]) + F[95])),-((-1 + d)*F[96]),OneMinusdd*F[97],d*(F[66] + F[96]) + F[98],-((-1 + d)*(d*(F[67] + F[97]) + F[99])),
        -((-1 + d)*F[100]),OneMinusdd*F[101],d*(F[70] + F[100]) + F[102],-((-1 + d)*(d*(F[71] + F[101]) + F[103])),OneMinusdd*F[104],OneMinusddd*F[105],-((-1 + d)*(d*(F[74] + F[104]) + F[106])),
        OneMinusdd*(d*(F[75] + F[105]) + F[107]),-((-1 + d)*(d*(F[44] + F[104]) + F[108])),OneMinusdd*(d*(F[45] + F[105]) + F[109]),
        dd*(F[14] + F[44] + F[74] + F[104]) + d*(F[46] + F[78] + F[106] + F[108]) + F[110],-((-1 + d)*(dd*(F[15] + F[45] + F[75] + F[105]) + d*(F[47] + F[79] + F[107] + F[109]) + F[111])),OneMinusdd*F[112],
        -((-1 + d)*(d*(F[97] + F[112]) + F[113])),-((-1 + d)*(d*(F[82] + F[112]) + F[114])),dd*(F[67] + F[82] + F[97] + F[112]) + d*(F[83] + F[99] + F[113] + F[114]) + F[115],OneMinusdd*F[116],
        -((-1 + d)*(d*(F[101] + F[116]) + F[117])),-((-1 + d)*(d*(F[86] + F[116]) + F[118])),dd*(F[71] + F[86] + F[101] + F[116]) + d*(F[87] + F[103] + F[117] + F[118]) + F[119],OneMinusddd*F[120],
        OneMinusdd*(d*(F[105] + F[120]) + F[121]),OneMinusdd*(d*(F[90] + F[120]) + F[122]),-((-1 + d)*(dd*(F[75] + F[90] + F[105] + F[120]) + d*(F[91] + F[107] + F[121] + F[122]) + F[123])),
        OneMinusdd*(d*(F[60] + F[120]) + F[124]),-((-1 + d)*(dd*(F[45] + F[60] + F[105] + F[120]) + d*(F[61] + F[109] + F[121] + F[124]) + F[125])),
        -((-1 + d)*(dd*(F[30] + F[60] + F[90] + F[120]) + d*(F[62] + F[94] + F[122] + F[124]) + F[126])),
        d*F[63] + d*F[95] + d*F[111] + ddd*(F[15] + F[30] + F[45] + F[60] + F[75] + F[90] + F[105] + F[120]) + d*F[123] + 
            dd*(F[31] + F[47] + F[61] + F[62] + F[79] + F[91] + F[94] + F[107] + F[109] + F[121] + F[122] + F[124]) + d*F[125] + d*F[126] + F[127],F[128],F[129],F[130],F[131],-((-1 + d)*F[132]),
        OneMinusdd*F[133],OneMinusdd*F[134],OneMinusddd*F[135],F[136],F[137],F[138],F[139],-((-1 + d)*F[140]),OneMinusdd*F[141],OneMinusdd*F[142],OneMinusddd*F[143],F[144],F[145],F[146],F[147],
        OneMinusdd*F[148],-((-1 + d)*(d*(F[133] + F[148]) + F[149])),OneMinusddd*F[150],OneMinusdd*(d*(F[135] + F[150]) + F[151]),F[152],F[153],F[154],F[155],OneMinusdd*F[156],
        -((-1 + d)*(d*(F[141] + F[156]) + F[157])),OneMinusddd*F[158],OneMinusdd*(d*(F[143] + F[158]) + F[159]),F[160],F[161],F[162],F[163],OneMinusdd*F[164],OneMinusddd*F[165],
        -((-1 + d)*(d*(F[134] + F[164]) + F[166])),OneMinusdd*(d*(F[135] + F[165]) + F[167]),F[168],F[169],F[170],F[171],OneMinusdd*F[172],OneMinusddd*F[173],-((-1 + d)*(d*(F[142] + F[172]) + F[174])),
        OneMinusdd*(d*(F[143] + F[173]) + F[175]),F[176],F[177],F[178],F[179],OneMinusddd*F[180],OneMinusdd*(d*(F[165] + F[180]) + F[181]),OneMinusdd*(d*(F[150] + F[180]) + F[182]),
        -((-1 + d)*(dd*(F[135] + F[150] + F[165] + F[180]) + d*(F[151] + F[167] + F[181] + F[182]) + F[183])),F[184],F[185],F[186],F[187],OneMinusddd*F[188],OneMinusdd*(d*(F[173] + F[188]) + F[189]),
        OneMinusdd*(d*(F[158] + F[188]) + F[190]),-((-1 + d)*(dd*(F[143] + F[158] + F[173] + F[188]) + d*(F[159] + F[175] + F[189] + F[190]) + F[191])),-((-1 + d)*F[192]),OneMinusdd*F[193],OneMinusdd*F[194],
        OneMinusddd*F[195],d*(F[132] + F[192]) + F[196],-((-1 + d)*(d*(F[133] + F[193]) + F[197])),-((-1 + d)*(d*(F[134] + F[194]) + F[198])),OneMinusdd*(d*(F[135] + F[195]) + F[199]),-((-1 + d)*F[200]),
        OneMinusdd*F[201],OneMinusdd*F[202],OneMinusddd*F[203],d*(F[140] + F[200]) + F[204],-((-1 + d)*(d*(F[141] + F[201]) + F[205])),-((-1 + d)*(d*(F[142] + F[202]) + F[206])),
        OneMinusdd*(d*(F[143] + F[203]) + F[207]),OneMinusdd*F[208],-((-1 + d)*(d*(F[193] + F[208]) + F[209])),OneMinusddd*F[210],OneMinusdd*(d*(F[195] + F[210]) + F[211]),
        -((-1 + d)*(d*(F[148] + F[208]) + F[212])),dd*(F[133] + F[148] + F[193] + F[208]) + d*(F[149] + F[197] + F[209] + F[212]) + F[213],OneMinusdd*(d*(F[150] + F[210]) + F[214]),
        -((-1 + d)*(dd*(F[135] + F[150] + F[195] + F[210]) + d*(F[151] + F[199] + F[211] + F[214]) + F[215])),OneMinusdd*F[216],-((-1 + d)*(d*(F[201] + F[216]) + F[217])),OneMinusddd*F[218],
        OneMinusdd*(d*(F[203] + F[218]) + F[219]),-((-1 + d)*(d*(F[156] + F[216]) + F[220])),dd*(F[141] + F[156] + F[201] + F[216]) + d*(F[157] + F[205] + F[217] + F[220]) + F[221],
        OneMinusdd*(d*(F[158] + F[218]) + F[222]),-((-1 + d)*(dd*(F[143] + F[158] + F[203] + F[218]) + d*(F[159] + F[207] + F[219] + F[222]) + F[223])),OneMinusdd*F[224],OneMinusddd*F[225],
        -((-1 + d)*(d*(F[194] + F[224]) + F[226])),OneMinusdd*(d*(F[195] + F[225]) + F[227]),-((-1 + d)*(d*(F[164] + F[224]) + F[228])),OneMinusdd*(d*(F[165] + F[225]) + F[229]),
        dd*(F[134] + F[164] + F[194] + F[224]) + d*(F[166] + F[198] + F[226] + F[228]) + F[230],-((-1 + d)*(dd*(F[135] + F[165] + F[195] + F[225]) + d*(F[167] + F[199] + F[227] + F[229]) + F[231])),
        OneMinusdd*F[232],OneMinusddd*F[233],-((-1 + d)*(d*(F[202] + F[232]) + F[234])),OneMinusdd*(d*(F[203] + F[233]) + F[235]),-((-1 + d)*(d*(F[172] + F[232]) + F[236])),
        OneMinusdd*(d*(F[173] + F[233]) + F[237]),dd*(F[142] + F[172] + F[202] + F[232]) + d*(F[174] + F[206] + F[234] + F[236]) + F[238],
        -((-1 + d)*(dd*(F[143] + F[173] + F[203] + F[233]) + d*(F[175] + F[207] + F[235] + F[237]) + F[239])),OneMinusddd*F[240],OneMinusdd*(d*(F[225] + F[240]) + F[241]),
        OneMinusdd*(d*(F[210] + F[240]) + F[242]),-((-1 + d)*(dd*(F[195] + F[210] + F[225] + F[240]) + d*(F[211] + F[227] + F[241] + F[242]) + F[243])),OneMinusdd*(d*(F[180] + F[240]) + F[244]),
        -((-1 + d)*(dd*(F[165] + F[180] + F[225] + F[240]) + d*(F[181] + F[229] + F[241] + F[244]) + F[245])),
        -((-1 + d)*(dd*(F[150] + F[180] + F[210] + F[240]) + d*(F[182] + F[214] + F[242] + F[244]) + F[246])),
        d*F[183] + d*F[215] + d*F[231] + ddd*(F[135] + F[150] + F[165] + F[180] + F[195] + F[210] + F[225] + F[240]) + d*F[243] + 
            dd*(F[151] + F[167] + F[181] + F[182] + F[199] + F[211] + F[214] + F[227] + F[229] + F[241] + F[242] + F[244]) + d*F[245] + d*F[246] + F[247],OneMinusddd*F[248],
        OneMinusdd*(d*(F[233] + F[248]) + F[249]),OneMinusdd*(d*(F[218] + F[248]) + F[250]),-((-1 + d)*(dd*(F[203] + F[218] + F[233] + F[248]) + d*(F[219] + F[235] + F[249] + F[250]) + F[251])),
        OneMinusdd*(d*(F[188] + F[248]) + F[252]),-((-1 + d)*(dd*(F[173] + F[188] + F[233] + F[248]) + d*(F[189] + F[237] + F[249] + F[252]) + F[253])),
        -((-1 + d)*(dd*(F[158] + F[188] + F[218] + F[248]) + d*(F[190] + F[222] + F[250] + F[252]) + F[254])),
        d*F[191] + d*F[223] + d*F[239] + ddd*(F[143] + F[158] + F[173] + F[188] + F[203] + F[218] + F[233] + F[248]) + d*F[251] + 
            dd*(F[159] + F[175] + F[189] + F[190] + F[207] + F[219] + F[222] + F[235] + F[237] + F[249] + F[250] + F[252]) + d*F[253] + d*F[254] + F[255] 
    };  

    std::vector<double> F31 = recomb1(F2,r);
    std::vector<double> F32 = recomb2(F2,r);
    std::vector<double> F33 = recomb3(F2,r);
    std::vector<double> F34 = recomb4(F2,r);

    F31.insert( F31.end(), F32.begin(), F32.end() );
    F31.insert( F31.end(), F33.begin(), F33.end() );
    F31.insert( F31.end(), F34.begin(), F34.end() );

    return F31; 
}

std::vector<int> initialize_pop(const int &n_wild, const int &n_drive){
    std::vector<int> pop(n_genotypes,0);
    pop[0] = n_wild;
    pop[n_genotypes-1] = n_drive;
    return pop;
}

rnd::discrete_distribution make_dist(const std::vector<int> &pop){
    rnd::discrete_distribution dist(n_genotypes);
    for(int i = 0; i < n_genotypes;++i){
        dist[i] = pop[i];
    }
    return dist;
}

const std::vector<double> initialize_multiplicative_payload(const double &sd, const double &spl, const double &e){
    const std::vector<double> payload{
                1,(1 - e)*(1 - spl),(1 - e)*(1 - spl),1 - spl,1 - sd,(1 - e)*(1 - \
        sd)*(1 - spl),(1 - e)*(1 - sd)*(1 - spl),(1 - sd)*(1 - spl),1 - sd,(1 \
        - e)*(1 - sd)*(1 - spl),(1 - e)*(1 - sd)*(1 - spl),(1 - sd)*(1 - \
        spl),pow(1 - sd,2),(1 - e)*pow(1 - sd,2)*(1 - spl),(1 - e)*pow(1 - \
        sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),(1 - e)*(1 - spl),(1 - e)*(1 \
        - spl),1 - spl,1 - spl,(1 - e)*(1 - sd)*(1 - spl),(1 - e)*(1 - sd)*(1 \
        - spl),(1 - sd)*(1 - spl),(1 - sd)*(1 - spl),(1 - e)*(1 - sd)*(1 - \
        spl),(1 - e)*(1 - sd)*(1 - spl),(1 - sd)*(1 - spl),(1 - sd)*(1 - \
        spl),(1 - e)*pow(1 - sd,2)*(1 - spl),(1 - e)*pow(1 - sd,2)*(1 - \
        spl),pow(1 - sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),(1 - e)*(1 - \
        spl),1 - spl,(1 - e)*(1 - spl),1 - spl,(1 - e)*(1 - sd)*(1 - spl),(1 \
        - sd)*(1 - spl),(1 - e)*(1 - sd)*(1 - spl),(1 - sd)*(1 - spl),(1 - \
        e)*(1 - sd)*(1 - spl),(1 - sd)*(1 - spl),(1 - e)*(1 - sd)*(1 - \
        spl),(1 - sd)*(1 - spl),(1 - e)*pow(1 - sd,2)*(1 - spl),pow(1 - \
        sd,2)*(1 - spl),(1 - e)*pow(1 - sd,2)*(1 - spl),pow(1 - sd,2)*(1 - \
        spl),1 - spl,1 - spl,1 - spl,1 - spl,(1 - sd)*(1 - spl),(1 - sd)*(1 - \
        spl),(1 - sd)*(1 - spl),(1 - sd)*(1 - spl),(1 - sd)*(1 - spl),(1 - \
        sd)*(1 - spl),(1 - sd)*(1 - spl),(1 - sd)*(1 - spl),pow(1 - sd,2)*(1 \
        - spl),pow(1 - sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),pow(1 - \
        sd,2)*(1 - spl),1 - sd,(1 - e)*(1 - sd)*(1 - spl),(1 - e)*(1 - sd)*(1 \
        - spl),(1 - sd)*(1 - spl),pow(1 - sd,2),(1 - e)*pow(1 - sd,2)*(1 - \
        spl),(1 - e)*pow(1 - sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),pow(1 - \
        sd,2),(1 - e)*pow(1 - sd,2)*(1 - spl),(1 - e)*pow(1 - sd,2)*(1 - \
        spl),pow(1 - sd,2)*(1 - spl),pow(1 - sd,3),(1 - e)*pow(1 - sd,3)*(1 - \
        spl),(1 - e)*pow(1 - sd,3)*(1 - spl),pow(1 - sd,3)*(1 - spl),(1 - \
        e)*(1 - sd)*(1 - spl),(1 - e)*(1 - sd)*(1 - spl),(1 - sd)*(1 - \
        spl),(1 - sd)*(1 - spl),(1 - e)*pow(1 - sd,2)*(1 - spl),(1 - e)*pow(1 \
        - sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),(1 \
        - e)*pow(1 - sd,2)*(1 - spl),(1 - e)*pow(1 - sd,2)*(1 - spl),pow(1 - \
        sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),(1 - e)*pow(1 - sd,3)*(1 - \
        spl),(1 - e)*pow(1 - sd,3)*(1 - spl),pow(1 - sd,3)*(1 - spl),pow(1 - \
        sd,3)*(1 - spl),(1 - e)*(1 - sd)*(1 - spl),(1 - sd)*(1 - spl),(1 - \
        e)*(1 - sd)*(1 - spl),(1 - sd)*(1 - spl),(1 - e)*pow(1 - sd,2)*(1 - \
        spl),pow(1 - sd,2)*(1 - spl),(1 - e)*pow(1 - sd,2)*(1 - spl),pow(1 - \
        sd,2)*(1 - spl),(1 - e)*pow(1 - sd,2)*(1 - spl),pow(1 - sd,2)*(1 - \
        spl),(1 - e)*pow(1 - sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),(1 - \
        e)*pow(1 - sd,3)*(1 - spl),pow(1 - sd,3)*(1 - spl),(1 - e)*pow(1 - \
        sd,3)*(1 - spl),pow(1 - sd,3)*(1 - spl),(1 - sd)*(1 - spl),(1 - \
        sd)*(1 - spl),(1 - sd)*(1 - spl),(1 - sd)*(1 - spl),pow(1 - sd,2)*(1 \
        - spl),pow(1 - sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),pow(1 - \
        sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),pow(1 \
        - sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),pow(1 - sd,3)*(1 - \
        spl),pow(1 - sd,3)*(1 - spl),pow(1 - sd,3)*(1 - spl),pow(1 - sd,3)*(1 \
        - spl),1 - sd,(1 - e)*(1 - sd)*(1 - spl),(1 - e)*(1 - sd)*(1 - \
        spl),(1 - sd)*(1 - spl),pow(1 - sd,2),(1 - e)*pow(1 - sd,2)*(1 - \
        spl),(1 - e)*pow(1 - sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),pow(1 - \
        sd,2),(1 - e)*pow(1 - sd,2)*(1 - spl),(1 - e)*pow(1 - sd,2)*(1 - \
        spl),pow(1 - sd,2)*(1 - spl),pow(1 - sd,3),(1 - e)*pow(1 - sd,3)*(1 - \
        spl),(1 - e)*pow(1 - sd,3)*(1 - spl),pow(1 - sd,3)*(1 - spl),(1 - \
        e)*(1 - sd)*(1 - spl),(1 - e)*(1 - sd)*(1 - spl),(1 - sd)*(1 - \
        spl),(1 - sd)*(1 - spl),(1 - e)*pow(1 - sd,2)*(1 - spl),(1 - e)*pow(1 \
        - sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),(1 \
        - e)*pow(1 - sd,2)*(1 - spl),(1 - e)*pow(1 - sd,2)*(1 - spl),pow(1 - \
        sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),(1 - e)*pow(1 - sd,3)*(1 - \
        spl),(1 - e)*pow(1 - sd,3)*(1 - spl),pow(1 - sd,3)*(1 - spl),pow(1 - \
        sd,3)*(1 - spl),(1 - e)*(1 - sd)*(1 - spl),(1 - sd)*(1 - spl),(1 - \
        e)*(1 - sd)*(1 - spl),(1 - sd)*(1 - spl),(1 - e)*pow(1 - sd,2)*(1 - \
        spl),pow(1 - sd,2)*(1 - spl),(1 - e)*pow(1 - sd,2)*(1 - spl),pow(1 - \
        sd,2)*(1 - spl),(1 - e)*pow(1 - sd,2)*(1 - spl),pow(1 - sd,2)*(1 - \
        spl),(1 - e)*pow(1 - sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),(1 - \
        e)*pow(1 - sd,3)*(1 - spl),pow(1 - sd,3)*(1 - spl),(1 - e)*pow(1 - \
        sd,3)*(1 - spl),pow(1 - sd,3)*(1 - spl),(1 - sd)*(1 - spl),(1 - \
        sd)*(1 - spl),(1 - sd)*(1 - spl),(1 - sd)*(1 - spl),pow(1 - sd,2)*(1 \
        - spl),pow(1 - sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),pow(1 - \
        sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),pow(1 \
        - sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),pow(1 - sd,3)*(1 - \
        spl),pow(1 - sd,3)*(1 - spl),pow(1 - sd,3)*(1 - spl),pow(1 - sd,3)*(1 \
        - spl),pow(1 - sd,2),(1 - e)*pow(1 - sd,2)*(1 - spl),(1 - e)*pow(1 - \
        sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),pow(1 - sd,3),(1 - e)*pow(1 - \
        sd,3)*(1 - spl),(1 - e)*pow(1 - sd,3)*(1 - spl),pow(1 - sd,3)*(1 - \
        spl),pow(1 - sd,3),(1 - e)*pow(1 - sd,3)*(1 - spl),(1 - e)*pow(1 - \
        sd,3)*(1 - spl),pow(1 - sd,3)*(1 - spl),pow(1 - sd,4),(1 - e)*pow(1 - \
        sd,4)*(1 - spl),(1 - e)*pow(1 - sd,4)*(1 - spl),pow(1 - sd,4)*(1 - \
        spl),(1 - e)*pow(1 - sd,2)*(1 - spl),(1 - e)*pow(1 - sd,2)*(1 - \
        spl),pow(1 - sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),(1 - e)*pow(1 - \
        sd,3)*(1 - spl),(1 - e)*pow(1 - sd,3)*(1 - spl),pow(1 - sd,3)*(1 - \
        spl),pow(1 - sd,3)*(1 - spl),(1 - e)*pow(1 - sd,3)*(1 - spl),(1 - \
        e)*pow(1 - sd,3)*(1 - spl),pow(1 - sd,3)*(1 - spl),pow(1 - sd,3)*(1 - \
        spl),(1 - e)*pow(1 - sd,4)*(1 - spl),(1 - e)*pow(1 - sd,4)*(1 - \
        spl),pow(1 - sd,4)*(1 - spl),pow(1 - sd,4)*(1 - spl),(1 - e)*pow(1 - \
        sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),(1 - e)*pow(1 - sd,2)*(1 - \
        spl),pow(1 - sd,2)*(1 - spl),(1 - e)*pow(1 - sd,3)*(1 - spl),pow(1 - \
        sd,3)*(1 - spl),(1 - e)*pow(1 - sd,3)*(1 - spl),pow(1 - sd,3)*(1 - \
        spl),(1 - e)*pow(1 - sd,3)*(1 - spl),pow(1 - sd,3)*(1 - spl),(1 - \
        e)*pow(1 - sd,3)*(1 - spl),pow(1 - sd,3)*(1 - spl),(1 - e)*pow(1 - \
        sd,4)*(1 - spl),pow(1 - sd,4)*(1 - spl),(1 - e)*pow(1 - sd,4)*(1 - \
        spl),pow(1 - sd,4)*(1 - spl),pow(1 - sd,2)*(1 - spl),pow(1 - sd,2)*(1 \
        - spl),pow(1 - sd,2)*(1 - spl),pow(1 - sd,2)*(1 - spl),pow(1 - \
        sd,3)*(1 - spl),pow(1 - sd,3)*(1 - spl),pow(1 - sd,3)*(1 - spl),pow(1 \
        - sd,3)*(1 - spl),pow(1 - sd,3)*(1 - spl),pow(1 - sd,3)*(1 - \
        spl),pow(1 - sd,3)*(1 - spl),pow(1 - sd,3)*(1 - spl),pow(1 - sd,4)*(1 \
        - spl),pow(1 - sd,4)*(1 - spl),pow(1 - sd,4)*(1 - spl),pow(1 - \
        sd,4)*(1 - spl)
    };
    return payload;
}

void print(std::vector<int> &pop, std::ofstream &output_file){
    for(int i = 0; i < n_genotypes; ++i){
        output_file << pop[i] << output_file.fill();
    }
    output_file << '\n';
}

void print_header(std::ofstream &output_file){
    for(int i = 0; i < n_genotypes; ++i){
        output_file << "type" << i << output_file.fill();
    }
    output_file << '\n';
}