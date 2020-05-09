#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include "Ngram.h"
#include "VocabMap.h"
using namespace std;
static Vocab voc;
static Vocab vocZhuYin, vocBig5;

// Get P(W2 | W1) -- bigram
LogP getUnigramProb(const char *w, Ngram &lm)
{
    VocabIndex wid = voc.getIndex(w);
    if(wid == Vocab_None)
        wid  = voc.getIndex(Vocab_Unknown);
    VocabIndex context[] = {Vocab_None};
    return lm.wordProb(wid, context);
}
LogP getBigramProb(const char *w1, const char *w2, Ngram &lm)
{
    VocabIndex wid1 = voc.getIndex(w1);
    VocabIndex wid2 = voc.getIndex(w2);

    if(wid1 == Vocab_None)  //OOV
        wid1 = voc.getIndex(Vocab_Unknown);
    if(wid2 == Vocab_None)  //OOV
        wid2 = voc.getIndex(Vocab_Unknown);

    VocabIndex context[] = { wid1, Vocab_None };
    return lm.wordProb( wid2, context);
}

// Get P(W3 | W1, W2) -- trigram
LogP getTrigramProb(const char *w1, const char *w2, const char *w3, Ngram &lm) 
{
    VocabIndex wid1 = voc.getIndex(w1);
    VocabIndex wid2 = voc.getIndex(w2);
    VocabIndex wid3 = voc.getIndex(w3);

    if(wid1 == Vocab_None)  //OOV
        wid1 = voc.getIndex(Vocab_Unknown);
    if(wid2 == Vocab_None)  //OOV
        wid2 = voc.getIndex(Vocab_Unknown);
    if(wid3 == Vocab_None)  //OOV
        wid3 = voc.getIndex(Vocab_Unknown);

    VocabIndex context[] = { wid2, wid1, Vocab_None };
    return lm.wordProb( wid3, context );
}
int viterbi(VocabString words[], int numWords, VocabMap &map, Ngram &lm)
{  
    printf("<s> ");
    LogP probAll[numWords][1024];
    char candidates[numWords][1024][3];
    int traceBack[numWords][1024];
    int candidateSize[numWords];
    
    memset(candidates, 0, sizeof(candidates));
    memset(traceBack, 0, sizeof(traceBack));
    memset(candidateSize, 0, sizeof(candidateSize));
    for(int i = 0; i < numWords; ++i)
        for(int j = 0; j < 1024; ++j)
            probAll[i][j] = INT_MIN;
    {
        VocabMapIter initMapIter(map, vocZhuYin.getIndex(words[0]));
        VocabIndex wid;
        Prob p;
        for(int i = 0; initMapIter.next(wid, p); ++i)
        {
            strcpy(candidates[0][i], vocBig5.getWord(wid));
            probAll[0][i] = getBigramProb("<s>", candidates[0][i], lm);
            if(probAll[0][i] == LogP_Zero)
                probAll[0][i] = INT_MIN;
            ++candidateSize[0];
        }
    }
    for(int word_i = 1; word_i < numWords; ++word_i)
    {
        VocabMapIter ZhuYinMapIter(map, vocZhuYin.getIndex(words[word_i]));
        VocabIndex wid;
        Prob p;
        for(int j = 0; ZhuYinMapIter.next(wid, p); ++j)
        {
            ++candidateSize[word_i];
            strcpy(candidates[word_i][j], vocBig5.getWord(wid));
            for(int i = 0; i < candidateSize[word_i - 1]; ++i)
            {
                LogP p_bigram = getBigramProb( candidates[word_i - 1][i], candidates[word_i][j], lm );
                if(p_bigram == LogP_Zero)
                    p_bigram = INT_MIN;
                LogP p_current = p_bigram + probAll[word_i - 1][i];
                if(p_current > probAll[word_i][j]){
                    probAll[word_i][j] = p_current;
                    traceBack[word_i][j] = i;
                }
            }
        }
    }
    int traced_i[numWords] = {-1};
    LogP prob_max = INT_MIN;
    for(int i = 0; i < candidateSize[numWords - 1]; ++i)
    {
        LogP p_bigram = getBigramProb( candidates[numWords - 1][i], "</s>", lm );
        if(p_bigram == LogP_Zero)
            p_bigram = INT_MIN;
        LogP p_current = p_bigram + probAll[numWords - 1][i];
        if(p_current > prob_max){
            prob_max = p_current;
            traced_i[numWords - 1] = i;
        }
    }
    for(int i = numWords - 1; i > 0; --i)
        traced_i[i - 1] = traceBack[i][traced_i[i]];

    for(int i = 0; i < numWords; ++i)
        printf( "%s ", candidates[i][traced_i[i]] );
    printf("</s>\n");
}
int viterbi3(VocabString words[], const int numWords, VocabMap &map, Ngram &lm)
{  
    printf("<s> ");
    LogP probAll[numWords][1024];
    char candidates[numWords][1024][3];
    int traceBack[numWords][1024];
    int candidateSize[numWords];
    
    memset(candidates, 0, sizeof(candidates));
    memset(traceBack, 0, sizeof(traceBack));
    memset(candidateSize, 0, sizeof(candidateSize));
    for(int i = 0; i < numWords; ++i)
        for(int j = 0; j < 1024; ++j)
            probAll[i][j] = INT_MIN;
    {
        VocabMapIter initMapIter(map, vocZhuYin.getIndex(words[0]));
        VocabIndex wid;
        Prob p;
        for(int i = 0; initMapIter.next(wid, p); ++i)
        {
            strcpy(candidates[0][i], vocBig5.getWord(wid));
            probAll[0][i] = getTrigramProb(Vocab_Unknown, "<s>", candidates[0][i], lm);
            if(probAll[0][i] == LogP_Zero)
                probAll[0][i] = INT_MIN;
            ++candidateSize[0];
        }
    }
    {
        VocabMapIter ZhuYinMapIter(map, vocZhuYin.getIndex(words[1]));
        VocabIndex wid;
        Prob p;
        for(int j = 0; ZhuYinMapIter.next(wid, p); ++j)
        {
            ++candidateSize[1];
            strcpy(candidates[1][j], vocBig5.getWord(wid));
            for(int i = 0; i < candidateSize[0]; ++i)
            {
                LogP p_trigram = getTrigramProb( "<s>", candidates[0][i], candidates[1][j], lm );
                if(p_trigram == LogP_Zero)
                    p_trigram = INT_MIN;
                LogP p_current = p_trigram + probAll[0][i];
                if(p_current > probAll[1][j]){
                    probAll[1][j] = p_current;
                    traceBack[1][j] = i;
                }
            }
        }
    }
    for(int word_i = 2; word_i < numWords; ++word_i)
    {
        VocabMapIter ZhuYinMapIter(map, vocZhuYin.getIndex(words[word_i]));
        VocabIndex wid;
        Prob p;
        for(int k = 0; ZhuYinMapIter.next(wid, p); ++k)
        {
            ++candidateSize[word_i];
            strcpy(candidates[word_i][k], vocBig5.getWord(wid));
            for(int j = 0; j < candidateSize[word_i - 1]; ++j)
            {
                for(int i = 0; i < candidateSize[word_i - 2]; ++i)
                {
                    LogP p_trigram = getTrigramProb( candidates[word_i - 2][i], candidates[word_i - 1][j], candidates[word_i][k], lm );
                    if(p_trigram == LogP_Zero)
                        p_trigram = INT_MIN;
                    LogP p_current = p_trigram + probAll[word_i - 1][j];
                    if(p_current > probAll[word_i][k]){
                        probAll[word_i][k] = p_current;
                        traceBack[word_i][k] = j;
                    }
                }
            }
        }
    }
    int traced_i[numWords] = {-1};
    
    LogP prob_max = INT_MIN;
    for(int j = 0; j < candidateSize[numWords - 1]; ++j)
    {
        for(int i = 0; i < candidateSize[numWords - 2]; ++i)
        {
            LogP p_trigram = getTrigramProb( candidates[numWords - 2][i], candidates[numWords - 1][j], "</s>", lm );
            if(p_trigram == LogP_Zero)
                p_trigram = INT_MIN;
            LogP p_current = p_trigram + probAll[numWords - 1][j];
            if(p_current > prob_max){
                prob_max = p_current;
                traced_i[numWords - 1] = j;
            }
        }
    }
    for(int i = numWords - 1; i > 0; --i)
        traced_i[i - 1] = traceBack[i][traced_i[i]];

    for(int i = 0; i < numWords; ++i)
        printf( "%s ", candidates[i][traced_i[i]] );
    printf("</s>\n");
}
int main(int argc, char* argv[])
{
    char *file_name, *map_name, *LM_name;
    int ngram_order;
    assert(argc % 2);
    for(int i = 1; i < argc; ++i)
    {
        if(strcmp(argv[i], "-text") == 0) file_name = argv[++i];
        else if(strcmp(argv[i], "-map") == 0) map_name = argv[++i];
        else if(strcmp(argv[i], "-lm") == 0) LM_name = argv[++i];
        else if(strcmp(argv[i], "-order") == 0) ngram_order = atoi(argv[++i]);
    }
    Ngram lm(voc, ngram_order);
    VocabMap map(vocZhuYin, vocBig5);
    {
        File lmFile( LM_name, "r" );
        lm.read(lmFile);
        lmFile.close();

        File mapFile( map_name, "r");
        map.read(mapFile);
        mapFile.close();
    }
    File testFile( file_name, "r");
    char *line;
    while(line = testFile.getline())
    {
        VocabString words[maxWordsPerLine];
        int numWords = Vocab::parseWords(line, words, maxWordsPerLine);
        assert(numWords < maxWordsPerLine);
        switch(ngram_order)
        {
            case 2:
                viterbi(words, numWords, map, lm);
                break;
            case 3:
                viterbi3(words, numWords, map, lm);
                break;
            default:
                fprintf(stderr, "invalid ngram order\n");
        }

        
    }
    return 0;
}
    /*
    VocabIndex wid = voc.getIndex("Ê¨");
    if(wid == Vocab_None) {
        printf("No word with wid = %d\n", wid);
        printf("where Vocab_None is %d\n", Vocab_None);
    }

    wid = voc.getIndex("±wŽÍ");
    VocabIndex context[] = {voc.getIndex("Å}") , voc.getIndex("¬r"), Vocab_None};
    printf("Log Prob(±wŽÍ|¬r-Å}) = %f\n", lm.wordProb(wid, context));
    */