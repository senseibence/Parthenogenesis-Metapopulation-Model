100 rem Facultative Parthenogenesis Metapopulation model v1
110 rem parth 19 adds option to have no migration, i.e. no subpopulation
200 rem get run parameters from user
300 print "enter size of metapopulation"
400 input popsize
401 rem PopSize=100
500 print "enter maximum size of subpopulation"
600 input maxsubsize
601 rem MaxSubSize=50
700 print "enter reciprocal of migration rate [enter 0 for no migration; i.e. no subpopulation]"
800 input migrecip
801 rem migRecip=1000
810 if micrecip = 0 then migration = 0 else migration = 1/migrecip
900 rem Migration = 1/migRecip
1000 print "enter reciprocal of mutation rate"
1100 input mutrecip
1101 rem mutRecip=1000
1200 mutation = 1/mutrecip
1300 rem ? "enter recombination rate between Locus 1 and Locus 2"
1400 rem input Rec1
1401 rec1 = 0.5
1500 rem ? "enter recombination rate between Locus 2 and Locus 3”
1600 rem input Rec2
1601 rec2 = 0.5
1700 print "enter number of individuals encountered in main population"
1800 input numinds
1801 rem NumInds=4
1802 print "enter number of individuals encountered in sub population"
1803 input numindssub
1810 print "enter number of offspring per female"
1820 input maxrepro
1821 rem MaxRepro=10
1830 print "choose fitness dominance [1] or fitness overdominance [2]"
1840 input f
1841 rem f=1
1850 if f = 1 then overdom = 0 else overdom = 1
1900 print "enter number of generations in the run"
2000 input numgens
2001 rem NumGens=5
2100 print "enter name of output file"
2200 input outfile$
2201 rem Outfile$="outfile"
2210 rem ? "enter random seed"
2220 rem input RandomSeed
2230 randomseed = rnd(5000)
2260 randomize randomseed
2270 print "enter % reproduction of Parthenogenesis vs. sexual reproduction"
2280 input parthreduction
2290 parthrepro = int(parthreduction*maxrepro+0.5)
2300 print "Sexual reproduction produces ";maxrepro;" offspring. Parthenogenetic reproduction produces ";parthrepro;" offspring."
2310 print "Enter % reproduction for parthenogenetic-capable females reproducing sexually"
2320 input x
2325 parthpenalty = 1-x
2330 print "Hit any key followed by <return> to continue"
2380 input continue
2390 rem dimensioning
2400 dim sex(popsize*2)
2500 dim location(popsize)
2600 dim loc1allele1(popsize)
2700 dim loc1allele2(popsize)
2800 dim loc2allele1(popsize)
2900 dim loc2allele2(popsize)
3000 dim loc3allele1(popsize)
3100 dim loc3allele2(popsize)
3200 dim offspringalive(maxrepro*popsize)
3400 dim offspringsex(maxrepro*popsize)
3500 dim offspringlocation(maxrepro*popsize)
3600 dim offspringloc1allele1(maxrepro*popsize)
3700 dim offspringloc1allele2(maxrepro*popsize)
3800 dim offspringloc2allele1(maxrepro*popsize)
3900 dim offspringloc2allele2(maxrepro*popsize)
4000 dim offspringloc3allele1(maxrepro*popsize)
4100 dim offspringloc3allele2(maxrepro*popsize)
4200 dim popinds(popsize)
4300 dim subinds(popsize)
4400 dim meeting(popsize,numinds)
4500 dim sexmothers(popsize)
4600 dim sexfathers(popsize)
4700 dim parthmothers(popsize)
4710 dim score(maxrepro*popsize)
4720 dim originalindex(maxrepro*popsize)
4730 dim newoffspringalive(maxrepro*popsize)
4740 dim newoffspringsex(maxrepro*popsize)
4750 dim newoffspringlocation(maxrepro*popsize)
4760 dim newoffspringloc1allele1(maxrepro*popsize)
4770 dim newoffspringloc1allele2(maxrepro*popsize)
4780 dim newoffspringloc2allele1(maxrepro*popsize)
4790 dim newoffspringloc2allele2(maxrepro*popsize)
4800 dim newoffspringloc3allele1(maxrepro*popsize)
4810 dim newoffspringloc3allele2(maxrepro*popsize)
4815 dim loc1freq(2*numgens)
4820 dim loc2freq(2*numgens)
4830 dim loc3freq(2*numgens)
4840 dim mainpopcount(2*numgens)
4850 dim subpopcount(2*numgens)
4860 dim sexualoffspring(2*numgens)
4870 dim parthoffspring(2*numgens)
4880 dim sexuals(2*numgens)
4890 dim parths(2*numgens)
4900 rem keep track of the maximum value of locus 2 allele that has appeared
4910 maxparth = 0
5000 rem set initial arrays
5100 for a = 1 to popsize
5200 aa = rnd(1)
5300 if aa < 0.5 then sex(a) = 1 else sex(a) = 0
5400 location(a) = 0
5500 loc1allele1(a) = 1
5600 if overdom = 1 then loc1allele2(a) = 0 else loc1allele2(a) = 1
5700 loc2allele1(a) = 0
5800 loc2allele2(a) = 0
5900 loc3allele1(a) = 0
6000 loc3allele2(a) = 0
6100 next a
6200 rem start main loop
6300 for main = 1 to numgens
6400 rem mutation phase
6500 for a = 1 to popsize
6600 rem locus 1 allele 1
6700 z = rnd(1)
6800 if z < mutation then loc1allele1(a) = abs(loc1allele1(a)-1)
6900 rem locus 1 allele 2
7000 z = rnd(1)
7100 if z < mutation then loc1allele2(a) = abs(loc1allele2(a)-1)
7200 rem locus 2 allele 1
7300 z = rnd(1)
7400 if z < mutation then loc2allele1(a) = rnd(1)/2
7410 if loc2allele1(a) > maxparth then maxparth = loc2allele1(a)
7500 rem locus 2 allele 2
7600 z = rnd(1)
7700 if z < mutation then loc2allele2(a) = rnd(1)/2
7710 if loc2allele2(a) > maxparth then maxparth = loc2allele2(a)
7800 rem locus 3 allele 1
7900 z = rnd(1)
8000 if z < mutation then loc3allele1(a) = rnd(1)/2
8100 rem locus 3 allele 2
8200 if z < mutation then loc3allele2(a) = rnd(1)/2
8300 next a
10000 rem encounter phase
10010 rem each individual encounters NumInds or NumIndsSub individuals depending on its environment
10020 rem for programming simplicity, self-encounters and duplicate encounters count (bad luck)
10030 rem values of the Meeting variable are the subscripts of the encountered individual
10100 for a = 1 to popsize
10150 if location(a) = 0 then x = numinds else x = numindssub
10200 for b = 1 to x
10300 zz = rnd(popsize+1)
10500 if location(a) = location(zz) then goto 10600 else goto 10300
10600 meeting(a,b) = zz
10700 next b
10800 next a
11000 rem encounters check
11100 rem for a = 1 to PopSize
11200 rem if Sex(a)=1 then goto 11300 else goto 11800
11300 rem ? "female ";a; "meets ";
11400 rem for b = 1 to NumInds
11500 rem if Sex(Meeting(a,b))=1 then ? "female ";Meeting(a,b);" ";
11600 rem if Sex(Meeting(a,b))=0 then ? "male ";Meeting(a,b);" ";
11700 rem next b
11750 rem ?
11800 rem next a
13900 sexualmatingscount = 1
13950 parthmatingscount = 0
14000 rem find male that each female mates with
14050 rem for computational simplicity, each female will mate with the last male she encountered
14100 for a = 1 to popsize-1
14110 matingflag = 0
14150 if sex(a) = 1 then goto 14200 else goto 15100
14200 rem focal individual is female
14250 for b = 1 to numinds
14400 if sex(meeting(a,b)) = 0 then 
14450      matingflag = 1
14500      sexmothers(sexualmatingscount) = a
14600      sexfathers(sexualmatingscount) = meeting(a,b)
14650 endif
15000 next b
15010 if matingflag = 1 then goto 15070 else goto 15050
15050 parthmatingscount = parthmatingscount+1
15060 parthmothers(parthmatingscount) = a
15065 goto 15100
15070 sexualmatingscount = sexualmatingscount+1
15100 next a
15110 rem record number of sexual and parthenogenetic reproducing females this generation
15120 sexuals(main) = sexualmatingscount-1
15130 parths(main) = parthmatingscount
15200 rem check matings
15300 rem for a = 1 to SexualMatingsCount-1
15400 rem ? "female ";SexMothers(a); " mates with male ";SexFathers(a)
15500 rem next a
15600 rem for a = 1 to ParthMatingsCount
15700 rem ? "female ";ParthMothers(a); " attempts parthenogenesis"
15800 rem next a
16000 rem reproduction phase
16100 rem sexual reproduction
16200 offspringcount = 0
16300 for a = 1 to sexualmatingscount-1
16400 if loc2allele1(sexmothers(a))+loc2allele2(sexmothers(a)) = 0 then max = maxrepro else max = int(maxrepro*(1-parthpenalty)+0.5)
16410 for b = 1 to max
16450 offspringcount = offspringcount+1
16500 rem determine maternal haplotype
16600 x = rnd(2)
16700 y = rnd(1)
16800 z = rnd(1)
16900 if x = 0 and y > rec1 and z > rec2 then mhaplotype = 1
17000 if x = 0 and y > rec1 and z < rec2 then mhaplotype = 2
17100 if x = 0 and y < rec1 and z > rec2 then mhaplotype = 3
17200 if x = 0 and y < rec1 and z < rec2 then mhaplotype = 4
17300 if x = 1 and y > rec1 and z > rec2 then mhaplotype = 5
17400 if x = 1 and y > rec1 and z < rec2 then mhaplotype = 6
17500 if x = 1 and y < rec1 and z > rec2 then mhaplotype = 7
17600 if x = 1 and y < rec1 and z < rec2 then mhaplotype = 8
17700 if mhaplotype = 1 then 
17710     offspringloc1allele1(offspringcount) = loc1allele1(sexmothers(a))
17720     offspringloc2allele1(offspringcount) = loc2allele1(sexmothers(a))
17730     offspringloc3allele1(offspringcount) = loc3allele1(sexmothers(a))
17740 endif
18000 if mhaplotype = 2 then 
18010     offspringloc1allele1(offspringcount) = loc1allele1(sexmothers(a))
18020     offspringloc2allele1(offspringcount) = loc2allele1(sexmothers(a))
18030     offspringloc3allele1(offspringcount) = loc3allele2(sexmothers(a))
18040 endif
18300 if mhaplotype = 3 then 
18310     offspringloc1allele1(offspringcount) = loc1allele1(sexmothers(a))
18320     offspringloc2allele1(offspringcount) = loc2allele2(sexmothers(a))
18330     offspringloc3allele1(offspringcount) = loc3allele2(sexmothers(a))
18340 endif
18600 if mhaplotype = 4 then 
18610     offspringloc1allele1(offspringcount) = loc1allele1(sexmothers(a))
18620     offspringloc2allele1(offspringcount) = loc2allele2(sexmothers(a))
18630     offspringloc3allele1(offspringcount) = loc3allele1(sexmothers(a))
18640 endif
18900 if mhaplotype = 5 then 
18910     offspringloc1allele1(offspringcount) = loc1allele2(sexmothers(a))
18920     offspringloc2allele1(offspringcount) = loc2allele2(sexmothers(a))
18930     offspringloc3allele1(offspringcount) = loc3allele2(sexmothers(a))
18940 endif
19200 if mhaplotype = 6 then 
19210     offspringloc1allele1(offspringcount) = loc1allele2(sexmothers(a))
19220     offspringloc2allele1(offspringcount) = loc2allele2(sexmothers(a))
19230     offspringloc3allele1(offspringcount) = loc3allele1(sexmothers(a))
19240 endif
19500 if mhaplotype = 7 then 
19510     offspringloc1allele1(offspringcount) = loc1allele2(sexmothers(a))
19520     offspringloc2allele1(offspringcount) = loc2allele1(sexmothers(a))
19530     offspringloc3allele1(offspringcount) = loc3allele1(sexmothers(a))
19540 endif
19800 if mhaplotype = 8 then 
19810     offspringloc1allele1(offspringcount) = loc1allele2(sexmothers(a))
19820     offspringloc2allele1(offspringcount) = loc2allele1(sexmothers(a))
19830     offspringloc3allele1(offspringcount) = loc3allele2(sexmothers(a))
19840 endif
20100 rem determine paternal haplotype
20200 x = rnd(2)
20300 y = rnd(1)
20400 z = rnd(1)
20500 if x = 0 and y > rec1 and z > rec2 then phaplotype = 1
20600 if x = 0 and y > rec1 and z < rec2 then phaplotype = 2
20700 if x = 0 and y < rec1 and z > rec2 then phaplotype = 3
20800 if x = 0 and y < rec1 and z < rec2 then phaplotype = 4
20900 if x = 1 and y > rec1 and z > rec2 then phaplotype = 5
21000 if x = 1 and y > rec1 and z < rec2 then phaplotype = 6
21100 if x = 1 and y < rec1 and z > rec2 then phaplotype = 7
21200 if x = 1 and y < rec1 and z < rec2 then phaplotype = 8
21300 if phaplotype = 1 then 
21310     offspringloc1allele2(offspringcount) = loc1allele1(sexfathers(a))
21320     offspringloc2allele2(offspringcount) = loc2allele1(sexfathers(a))
21330     offspringloc3allele2(offspringcount) = loc3allele1(sexfathers(a))
21340 endif
21600 if phaplotype = 2 then 
21610     offspringloc1allele2(offspringcount) = loc1allele1(sexfathers(a))
21620     offspringloc2allele2(offspringcount) = loc2allele1(sexfathers(a))
21630     offspringloc3allele2(offspringcount) = loc3allele2(sexfathers(a))
21640 endif
21900 if phaplotype = 3 then 
21910     offspringloc1allele2(offspringcount) = loc1allele1(sexfathers(a))
21920     offspringloc2allele2(offspringcount) = loc2allele2(sexfathers(a))
21930     offspringloc3allele2(offspringcount) = loc3allele2(sexfathers(a))
21940 endif
22200 if phaplotype = 4 then 
22210     offspringloc1allele2(offspringcount) = loc1allele1(sexfathers(a))
22220     offspringloc2allele2(offspringcount) = loc2allele2(sexfathers(a))
22230     offspringloc3allele2(offspringcount) = loc3allele1(sexfathers(a))
22240 endif
22500 if phaplotype = 5 then 
22510     offspringloc1allele2(offspringcount) = loc1allele2(sexfathers(a))
22520     offspringloc2allele2(offspringcount) = loc2allele2(sexfathers(a))
22530     offspringloc3allele2(offspringcount) = loc3allele2(sexfathers(a))
22540 endif
22800 if phaplotype = 6 then 
22810     offspringloc1allele2(offspringcount) = loc1allele2(sexfathers(a))
22820     offspringloc2allele2(offspringcount) = loc2allele2(sexfathers(a))
22830     offspringloc3allele2(offspringcount) = loc3allele1(sexfathers(a))
22840 endif
23100 if phaplotype = 7 then 
23110     offspringloc1allele2(offspringcount) = loc1allele2(sexfathers(a))
23120     offspringloc2allele2(offspringcount) = loc2allele1(sexfathers(a))
23130     offspringloc3allele2(offspringcount) = loc3allele1(sexfathers(a))
23140 endif
23400 if phaplotype = 8 then 
23410     offspringloc1allele2(offspringcount) = loc1allele2(sexfathers(a))
23420     offspringloc2allele2(offspringcount) = loc2allele1(sexfathers(a))
23430     offspringloc3allele2(offspringcount) = loc3allele2(sexfathers(a))
23440 endif
23700 rem set other offspring parameters
23800 offspringalive(offspringcount) = 1
23900 zz = rnd(1)
24000 if zz < 0.5 then offspringsex(offspringcount) = 1 else offspringsex(offspringcount) = 0
24100 if location(sexmothers(a)) = 1 then offspringlocation(offspringcount) = 1 else offspringlocation(offspringcount) = 0
24200 next b
24300 next a
24350 rem save current OffspringCount number as the number of offspring produced by sexual reproduction this generation
24360 sexualoffspring(main) = offspringcount
24400 rem parthenogenic reproduction
24500 rem produce parthogenetic offspring by gamete cell duplication
24510 for a = 1 to parthmatingscount
24600 zz = rnd(1)
24700 if zz < (loc2allele1(parthmothers(a))+loc2allele2(parthmothers(a))) then goto 25000 else goto 31900
25000 rem female will reproduce parthenogenetic offspring
25100 for b = 1 to parthrepro
25150 offspringcount = offspringcount+1
25200 rem determine maternal haplotype
25300 x = rnd(2)
25400 y = rnd(1)
25500 z = rnd(1)
25600 if x = 0 and y > rec1 and z > rec2 then mhaplotype = 1
25700 if x = 0 and y > rec1 and z < rec2 then mhaplotype = 2
25800 if x = 0 and y < rec1 and z > rec2 then mhaplotype = 3
25900 if x = 0 and y < rec1 and z < rec2 then mhaplotype = 4
26000 if x = 1 and y > rec1 and z > rec2 then mhaplotype = 5
26100 if x = 1 and y > rec1 and z < rec2 then mhaplotype = 6
26200 if x = 1 and y < rec1 and z > rec2 then mhaplotype = 7
26300 if x = 1 and y < rec1 and z < rec2 then mhaplotype = 8
26400 rem parthenogenetic offspring represent duplicated maternal meiotic haplotypes
26500 if mhaplotype = 1 then 
26510     offspringloc1allele1(offspringcount) = loc1allele1(parthmothers(a))
26520     offspringloc2allele1(offspringcount) = loc2allele1(parthmothers(a))
26530     offspringloc3allele1(offspringcount) = loc3allele1(parthmothers(a))
26540     offspringloc1allele2(offspringcount) = loc1allele1(parthmothers(a))
26550     offspringloc2allele2(offspringcount) = loc2allele1(parthmothers(a))
26560     offspringloc3allele2(offspringcount) = loc3allele1(parthmothers(a))
26570 endif
27100 if mhaplotype = 2 then 
27110     offspringloc1allele1(offspringcount) = loc1allele1(parthmothers(a))
27120     offspringloc2allele1(offspringcount) = loc2allele1(parthmothers(a))
27130     offspringloc3allele1(offspringcount) = loc3allele2(parthmothers(a))
27140     offspringloc1allele2(offspringcount) = loc1allele1(parthmothers(a))
27150     offspringloc2allele2(offspringcount) = loc2allele1(parthmothers(a))
27160     offspringloc3allele2(offspringcount) = loc3allele2(parthmothers(a))
27170 endif
27700 if mhaplotype = 3 then 
27710     offspringloc1allele1(offspringcount) = loc1allele1(parthmothers(a))
27720     offspringloc2allele1(offspringcount) = loc2allele2(parthmothers(a))
27730     offspringloc3allele1(offspringcount) = loc3allele2(parthmothers(a))
27740     offspringloc1allele2(offspringcount) = loc1allele1(parthmothers(a))
27750     offspringloc2allele2(offspringcount) = loc2allele2(parthmothers(a))
27760     offspringloc3allele2(offspringcount) = loc3allele2(parthmothers(a))
27770 endif
28300 if mhaplotype = 4 then 
28310     offspringloc1allele1(offspringcount) = loc1allele1(parthmothers(a))
28320     offspringloc2allele1(offspringcount) = loc2allele2(parthmothers(a))
28330     offspringloc3allele1(offspringcount) = loc3allele1(parthmothers(a))
28340     offspringloc1allele2(offspringcount) = loc1allele1(parthmothers(a))
28350     offspringloc2allele2(offspringcount) = loc2allele2(parthmothers(a))
28360     offspringloc3allele2(offspringcount) = loc3allele1(parthmothers(a))
28370 endif
28900 if mhaplotype = 5 then 
28910     offspringloc1allele1(offspringcount) = loc1allele2(parthmothers(a))
28920     offspringloc2allele1(offspringcount) = loc2allele2(parthmothers(a))
28930     offspringloc3allele1(offspringcount) = loc3allele2(parthmothers(a))
28940     offspringloc1allele2(offspringcount) = loc1allele2(parthmothers(a))
28950     offspringloc2allele2(offspringcount) = loc2allele2(parthmothers(a))
28960     offspringloc3allele2(offspringcount) = loc3allele2(parthmothers(a))
28970 endif
29500 if mhaplotype = 6 then 
29510     offspringloc1allele1(offspringcount) = loc1allele2(parthmothers(a))
29520     offspringloc2allele1(offspringcount) = loc2allele2(parthmothers(a))
29530     offspringloc3allele1(offspringcount) = loc3allele1(parthmothers(a))
29540     offspringloc1allele2(offspringcount) = loc1allele2(parthmothers(a))
29550     offspringloc2allele2(offspringcount) = loc2allele2(parthmothers(a))
29560     offspringloc3allele2(offspringcount) = loc3allele1(parthmothers(a))
29570 endif
30100 if mhaplotype = 7 then 
30110     offspringloc1allele1(offspringcount) = loc1allele2(parthmothers(a))
30120     offspringloc2allele1(offspringcount) = loc2allele1(parthmothers(a))
30130     offspringloc3allele1(offspringcount) = loc3allele1(parthmothers(a))
30140     offspringloc1allele2(offspringcount) = loc1allele2(parthmothers(a))
30150     offspringloc2allele2(offspringcount) = loc2allele1(parthmothers(a))
30160     offspringloc3allele2(offspringcount) = loc3allele1(parthmothers(a))
30170 endif
30200 if mhaplotype = 8 then 
30210     offspringloc1allele1(offspringcount) = loc1allele2(parthmothers(a))
30220     offspringloc2allele1(offspringcount) = loc2allele1(parthmothers(a))
30230     offspringloc3allele1(offspringcount) = loc3allele2(parthmothers(a))
30240     offspringloc1allele2(offspringcount) = loc1allele2(parthmothers(a))
30250     offspringloc2allele2(offspringcount) = loc2allele1(parthmothers(a))
30260     offspringloc3allele2(offspringcount) = loc3allele2(parthmothers(a))
30270 endif
31300 rem set other offspring parameters
31400 offspringalive(offspringcount) = 1
31500 zz = rnd(1)
31600 if zz < 0.5 then offspringsex(offspringcount) = 1 else offspringsex(offspringcount) = 0
31700 if location(parthmothers(a)) = 1 then offspringlocation(offspringcount) = 1 else offspringlocation(offspringcount) = 0
31800 next b
31900 next a
31950 rem set number of offspring produced in parthenogenetic loop as the number of parthenogenetic offspring produced this generation
31960 parthoffspring(main) = offspringcount-sexualoffspring(main)
41000 rem randomly sort offspring pool
41050 rem assign RND (0,1) to each member of OffspringCount
41100 for a = 1 to offspringcount
41200 score(a) = rnd(1)
41250 originalindex(a) = a
41300 next a
41400 rem sort OffspringCount subscripts by the random number and get new subscript for each member
41500 rem Bubble Sort Algorithm
41600 for i = 1 to offspringcount-1
41700   for j = 1 to offspringcount-i
41800     if score(j) > score(j+1) then 
41900       rem Swap the elements in score array
42000       temp = score(j)
42100       score(j) = score(j+1)
42200       score(j+1) = temp
42300       rem Swap the corresponding indices in originalIndex array
42400       tempindex = originalindex(j)
42500       originalindex(j) = originalindex(j+1)
42600       originalindex(j+1) = tempindex
42700     endif
42800   next j
42900 next i
43000 rem Print the sorted array and original indices
43100 rem PRINT "Sorted Scores and Original Indices:"
43200 rem FOR i = 1 TO OffspringCount
43300 rem   PRINT "Score: "; score(i); " Original Index: "; originalIndex(i)
43400 rem NEXT i
43500 rem set new reordered offspring array
43600 for a = 1 to offspringcount
43700 newoffspringalive(a) = offspringalive(originalindex(a))
43800 newoffspringsex(a) = offspringsex(originalindex(a))
43900 newoffspringlocation(a) = offspringlocation(originalindex(a))
44000 newoffspringloc1allele1(a) = offspringloc1allele1(originalindex(a))
44100 newoffspringloc1allele2(a) = offspringloc1allele2(originalindex(a))
44200 newoffspringloc2allele1(a) = offspringloc2allele1(originalindex(a))
44300 newoffspringloc2allele2(a) = offspringloc2allele2(originalindex(a))
44400 newoffspringloc3allele1(a) = offspringloc3allele1(originalindex(a))
44500 newoffspringloc3allele2(a) = offspringloc3allele2(originalindex(a))
44600 next a
45000 rem selection phase
45100 for a = 1 to offspringcount
45200 if overdom = 0 then goto 45300 else goto 46000
45300 rem fitness dominance
45400 if newoffspringloc1allele1(a)+newoffspringloc1allele2(a) = 0 then newoffspringalive(a) = 0
45500 goto 47000
46000 rem fitness overdominance
46100 if newoffspringloc1allele1(a)+newoffspringloc1allele2(a) = 0 then newoffspringalive(a) = 0
46200 if newoffspringloc1allele1(a)+newoffspringloc1allele2(a) = 2 then newoffspringalive(a) = 0
47000 next a
48000 rem simultaneously determine migration and send living offspring into adult pool for next generation
48100 rem need to keep count of number of individuals in the main and sub populations
48200 maincount = 0
48300 subcount = 0
48350 offspringcounter = 0
48400 for a = 1 to popsize
48450 offspringcounter = offspringcounter+1
48460 if newoffspringalive(offspringcounter) = 0 then goto 48470 else goto 48500
48470 offspringcounter = offspringcounter+1
48480 goto 48460
48500 z = rnd(1)
48600 if z < migration then newoffspringlocation(offspringcounter) = abs(newoffspringlocation(offspringcounter)-1)
48700 if newoffspringlocation(offspringcounter) = 0 then maincount = maincount+1 else subcount = subcount+1
48900 if subcount > maxsubsize then newoffspringlocation(offspringcounter) = 0
49000 rem Migration determination is complete; now the individual goes into the adult pool for the next generation
49100 sex(a) = newoffspringsex(offspringcounter)
49200 location(a) = newoffspringlocation(offspringcounter)
49300 loc1allele1(a) = newoffspringloc1allele1(offspringcounter)
49400 loc1allele2(a) = newoffspringloc1allele2(offspringcounter)
49500 loc2allele1(a) = newoffspringloc2allele1(offspringcounter)
49600 loc2allele2(a) = newoffspringloc2allele2(offspringcounter)
49700 loc3allele1(a) = newoffspringloc3allele1(offspringcounter)
49800 loc3allele2(a) = newoffspringloc3allele2(offspringcounter)
49900 next a
49910 rem check genotypes
49920 rem for a = 1 to PopSize
49930 rem ? "ind# ";a;" ";
49940 rem ? "sex ";Sex(a);" ";
49950 rem ? "loc ";Location(a);" ";
49960 rem ? "genotype ";Loc1Allele1(a);Loc1Allele2(a);" ";Loc2Allele1(a);Loc2Allele2(a);" ";Loc3Allele1(a);Loc3Allele2(a)
49965 rem next a
49970 rem ? "hit any key then [return] to continue"
49980 rem input x$
50000 rem compute and store generation data
50100 rem calculate allele frequencies at each locus for this generation
50110 cls
50112 maincount = 0
50114 subcount = 0
50120 print "generation ";main
50150 loc1count = 0
50200 loc2count = 0
50300 loc3count = 0
50400 for a = 1 to popsize
50450 loc1count = loc1count+loc1allele1(a)+loc1allele2(a)
50500 loc2count = loc2count+loc2allele1(a)+loc2allele2(a)
50600 loc3count = loc3count+loc3allele1(a)+loc3allele2(a)
50610 if location(a) = 0 then maincount = maincount+1
50620 if location(a) = 1 then subcount = subcount+1
50700 next a
50750 loc1freq(main) = loc1count/(2*popsize)
50760 print "mean fitness (determined by Locus 1) = ";loc1freq(main)
50800 loc2freq(main) = loc2count/popsize
50900 print "Mean Parthenogenetic Capability (determined by Locus 2) = ";loc2freq(main)
51000 loc3freq(main) = loc3count/popsize
51100 print "Locus 3 mean (neutral locus) = ";loc3freq(main)
51110 print "main population: ";maincount;" individuals"
51120 print "subpopulation: ";subcount;" individuals"
51125 print sexuals(main);" females reproduced sexually"
51130 print sexualoffspring(main);" sexual offspring produced"
51135 print parths(main);" females attempted to reproduce parthenogenetically"
51140 print parthoffspring(main);" parthenogenetic offspring produced ";
51142 if overdom = 1 and parthoffspring(main) > 0 then print "(but these die due to homoz. Lethality)" else print
51150 print "Maximum value of parthenogenetic allele that has appeared: ";maxparth
51200 next main
51300 rem print all data to output file
51400 open outfile$ for output as #1
51500 print #1,"metapopulation size: ";popsize
51600 print #1,"subpopulation size: ";maxsubsize
51700 print #1,"sexual reproduction produces ";maxrepro;" offspring"
51800 print #1,"parthenogenetic reproduction produces ";parthrepro;" offspring"
51850 print #1,"penalty for parthenogenetic-capable females reproducing sexually: ";parthpenalty
51900 print #1,"migration rate: ";migration
52000 print #1,"mutation rate: ";mutation
52100 print #1,"recombination rate, interval 1: ";rec1
52200 print #1,"recombination rate, interval 2: ";rec2
52300 print #1,"each individual in the main population encounters ";numinds;" individuals"
52400 print #1,"each individual in the subpopulation encounters ";numindssub;" individuals"
52500 print #1,"fitness overdominance (0=no, 1=yes): ";overdom
52510 print #1,"maximum value of parthenogenetic allele that appeared: ";maxparth
52550 print #1,"random seed: ";randomseed
52600 print #1,"generation";chr$(9);"mean Loc2";chr$(9);"mean Loc3";chr$(9);"Mean Loc1";chr$(9);"f prod sexually";chr$(9);"sexual offspring";chr$(9);"f prod parth";chr$(9);"parth offspring"
52700 for a = 1 to numgens
52800 print #1,a;chr$(9);loc2freq(a);chr$(9);loc3freq(a);chr$(9);loc1freq(a);chr$(9);sexuals(a);chr$(9);sexualoffspring(a);chr$(9);parths(a);chr$(9);parthoffspring(a)
52900 next a
53000 close #1
54000 end
