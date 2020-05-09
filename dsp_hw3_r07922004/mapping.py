# encoding=big5-hkscs
import sys
def main():
    if len(sys.argv) != 3:
        print("python3 mapping $(from) $(to)")
        return
    fpBig_ZhuYin = open(sys.argv[1], 'r', encoding='big5-hkscs')
    fpZhuYin_Big = open(sys.argv[2], 'w', encoding='big5-hkscs')
    if( (fpBig_ZhuYin is None) or (fpZhuYin_Big is None)):
        print("file not found or cannot create file")
        return
    zhuBigMap = {}
    for line in fpBig_ZhuYin.readlines():
        character, pronounce = line.split()
        zhuYins =  {zhuYin[0] for zhuYin in pronounce.split('/')} #set: noduplicated list
        for zhuYin in zhuYins:
            if zhuYin in zhuBigMap:
                zhuBigMap[zhuYin].append(character)
            else:
                zhuBigMap[zhuYin] = [character]
            zhuBigMap[character] = [character]
    fpBig_ZhuYin.close()
    # write the map
    for zhuYin, characters in zhuBigMap.items():
        buffer = zhuYin + '\t'
        buffer += " ".join(characters)
        buffer += '\n'
        fpZhuYin_Big.write(buffer)
    
    return

if __name__ == "__main__":
    main()