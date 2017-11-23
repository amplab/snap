using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;

namespace ErrorSim
{
    class Program
    {
        static double[] underlyingASEDistribution =
            {
                0.061375694
                ,0.090576348
                ,0.124467787
                ,0.15925689
                ,0.198495789
                ,0.238356875
                ,0.274425388
                ,0.311711966
                ,0.338824056
                ,0.380572883
                ,0.404787498
                ,0.432614617
                ,0.461148952
                ,0.483660718
                ,0.512669888
                ,0.534948601
                ,0.561709952
                ,0.58204263
                ,0.600275391
                ,0.638934588
                ,0.646386849
                ,0.661593218
                ,0.675332415
                ,0.696144384
                ,0.704558645
                ,0.722727007
                ,0.734169309
                ,0.751836223
                ,0.766332018
                ,0.776577984
                ,0.788065394
                ,0.795711704
                ,0.799108795
                ,0.823317062
                ,0.827050182
                ,0.832736479
                ,0.841376956
                ,0.848842158
                ,0.857585305
                ,0.866831427
                ,0.875172562
                ,0.881278573
                ,0.889126255
                ,0.892103815
                ,0.896176931
                ,0.904516966
                ,0.909361574
                ,0.913023276
                ,0.91534581
                ,0.91605211
                ,0.925432766
                ,0.927031296
                ,0.930831256
                ,0.934979879
                ,0.93665001
                ,0.938773245
                ,0.94021362
                ,0.94386287
                ,0.945187634
                ,0.945596911
                ,0.95304374
                ,0.954421059
                ,0.956191053
                ,0.959957013
                ,0.961346052
                ,0.962265629
                ,0.965476302
                ,0.966183152
                ,0.967031494
                ,0.968931931
                ,0.969856818
                ,0.971540195
                ,0.972253272
                ,0.973513088
                ,0.973847957
                ,0.97493863
                ,0.9758872
                ,0.976595943
                ,0.977207752
                ,0.97750679
                ,0.980262345
                ,0.982044974
                ,0.982477142
                ,0.983729694
                ,0.984795097
                ,0.985666451
                ,0.986368173
                ,0.986957032
                ,0.987644654
                ,0.988132857
                ,0.988691622
                ,0.9891291
                ,0.989549853
                ,0.989939841
                ,0.990257801
                ,0.990570451
                ,0.990839762
                ,0.991075135
                ,0.991264116
                ,0.991395598
                ,1
            };

        static double[] tumorGermlineRNAReadDepthDistribution =
        {
            0
            ,0
            ,0
            ,0
            ,0
            ,0
            ,0
            ,0
            ,0
            ,0
            ,0.032228023
            ,0.061294815
            ,0.08778307
            ,0.112046517
            ,0.134402081
            ,0.155199583
            ,0.174542306
            ,0.192737404
            ,0.209769191
            ,0.225861136
            ,0.242873024
            ,0.259061596
            ,0.274495219
            ,0.289267772
            ,0.303394454
            ,0.31689718
            ,0.329833938
            ,0.342270956
            ,0.354236558
            ,0.36573642
            ,0.376845806
            ,0.38754799
            ,0.397875568
            ,0.407892205
            ,0.41754327
            ,0.426841643
            ,0.435894816
            ,0.444640222
            ,0.453112105
            ,0.461661327
            ,0.469936354
            ,0.477994382
            ,0.485860316
            ,0.493448453
            ,0.500847854
            ,0.508034833
            ,0.515041681
            ,0.521854298
            ,0.528494659
            ,0.535006529
            ,0.541358473
            ,0.547494943
            ,0.553507319
            ,0.559406646
            ,0.565123707
            ,0.570724669
            ,0.576163995
            ,0.581504374
            ,0.586850247
            ,0.592029052
            ,0.597145047
            ,0.60213371
            ,0.607001086
            ,0.611775619
            ,0.616415741
            ,0.620972604
            ,0.625461711
            ,0.629841373
            ,0.634119402
            ,0.638336024
            ,0.642445642
            ,0.646461929
            ,0.65044849
            ,0.654341902
            ,0.658133316
            ,0.661869976
            ,0.665512939
            ,0.669179036
            ,0.672762118
            ,0.676273966
            ,0.679702981
            ,0.683078953
            ,0.686399866
            ,0.689645454
            ,0.692831712
            ,0.695950092
            ,0.699027819
            ,0.702044566
            ,0.705004485
            ,0.707900556
            ,0.710769829
            ,0.713599365
            ,0.716372378
            ,0.719096009
            ,0.721779658
            ,0.724426743
            ,0.727051488
            ,0.729638999
            ,0.732182499
            ,0.734674236
            ,0.73712648
            ,0.739529343
            ,0.741896436
            ,0.744224279
            ,0.746525936
            ,0.748788344
            ,0.751007536
            ,0.7532216
            ,0.755382864
            ,0.757520016
            ,0.759627259
            ,0.761693422
            ,0.763730773
            ,0.76574029
            ,0.767735036
            ,0.769693767
            ,0.77164145
            ,0.773556537
            ,0.775449649
            ,0.777309739
            ,0.779141323
            ,0.780942813
            ,0.782720864
            ,0.784463024
            ,0.786187969
            ,0.787889414
            ,0.789578102
            ,0.791244388
            ,0.792890652
            ,0.794507068
            ,0.796092719
            ,0.797654015
            ,0.799207498
            ,0.80075097
            ,0.802279609
            ,0.803771868
            ,0.805263822
            ,0.806718297
            ,0.808159832
            ,0.809583298
            ,0.810981432
            ,0.812361743
            ,0.813720079
            ,0.815071151
            ,0.816433028
            ,0.817746194
            ,0.819064182
            ,0.820352566
            ,0.821624956
            ,0.822889351
            ,0.824137265
            ,0.82535478
            ,0.826551298
            ,0.827753431
            ,0.828934322
            ,0.83011448
            ,0.831278219
            ,0.832428529
            ,0.833559428
            ,0.834672136
            ,0.835795832
            ,0.836893159
            ,0.837982489
            ,0.839063518
            ,0.84011952
            ,0.841170456
            ,0.842207047
            ,0.84323259
            ,0.844250563
            ,0.84524992
            ,0.846232612
            ,0.847222751
            ,0.848187986
            ,0.849151084
            ,0.850109788
            ,0.851038154
            ,0.851960599
            ,0.852891102
            ,0.853804879
            ,0.854696499
            ,0.855588302
            ,0.856473329
            ,0.85735213
            ,0.85822318
            ,0.859069446
            ,0.859910646
            ,0.860756913
            ,0.861588469
            ,0.862418926
            ,0.863236564
            ,0.864042361
            ,0.864843336
            ,0.865640404
            ,0.866418671
            ,0.867201516
            ,0.867971604
            ,0.868737786
            ,0.869486204
            ,0.870232486
            ,0.870973702
            ,0.871709852
            ,0.872433915
            ,0.873155964
            ,0.873865256
            ,0.874566124
            ,0.875250816
            ,0.875949791
            ,0.876647974
            ,0.877336572
            ,0.878003929
            ,0.878668356
            ,0.87932912
            ,0.879979568
            ,0.880630566
            ,0.881271798
            ,0.881912724
            ,0.88254663
            ,0.883175836
            ,0.883793689
            ,0.884412274
            ,0.885009679
            ,0.885618009
            ,0.886211995
            ,0.886789683
            ,0.887369446
            ,0.887959709
            ,0.888527508
            ,0.889094514
            ,0.889664877
            ,0.890230663
            ,0.890790527
            ,0.891344775
            ,0.891891882
            ,0.892427879
            ,0.892969859
            ,0.893500912
            ,0.894034406
            ,0.894575348
            ,0.895102372
            ,0.895614808
            ,0.896129441
            ,0.896637299
            ,0.897136672
            ,0.897628903
            ,0.898125773
            ,0.89861715
            ,0.899103338
            ,0.899589893
            ,0.900080354
            ,0.900552503
            ,0.901029108
            ,0.901492833
            ,0.90195723
            ,0.902409113
            ,0.902861607
            ,0.903323806
            ,0.903757194
            ,0.904215731
            ,0.904658031
            ,0.905100881
            ,0.905533963
            ,0.905969427
            ,0.906399153
            ,0.906819417
            ,0.907238705
            ,0.907647188
            ,0.908065072
            ,0.9084738
            ,0.908884481
            ,0.909292598
            ,0.909695649
            ,0.910098944
            ,0.910499065
            ,0.910889847
            ,0.911281178
            ,0.911662316
            ,0.912042476
            ,0.912425811
            ,0.912803469
            ,0.913173192
            ,0.913548653
            ,0.913916666
            ,0.914275402
            ,0.914646223
            ,0.915007584
            ,0.915358445
            ,0.915714678
            ,0.916068347
            ,0.916419758
            ,0.916766285
            ,0.917102863
            ,0.91745586
            ,0.917810323
            ,0.918154592
            ,0.918485249
            ,0.918810229
            ,0.91913881
            ,0.919474411
            ,0.91980635
            ,0.920129194
            ,0.920446544
            ,0.92076426
            ,0.921077703
            ,0.92139896
            ,0.921718385
            ,0.922032378
            ,0.922338741
            ,0.922647484
            ,0.922950123
            ,0.923252091
            ,0.923557111
            ,0.923854622
            ,0.924154332
            ,0.924447204
            ,0.924738246
            ,0.925028677
            ,0.925315079
            ,0.925604595
            ,0.92589588
            ,0.926176423
            ,0.926465999
            ,0.926742391
            ,0.927022873
            ,0.927297372
            ,0.927565462
            ,0.927838009
            ,0.928108785
            ,0.928376142
            ,0.92864936
            ,0.92890982
            ,0.929174919
            ,0.929442155
            ,0.929704874
            ,0.929957215
            ,0.930213525
            ,0.930464096
            ,0.930722481
            ,0.930972564
            ,0.931228324
            ,0.931476943
            ,0.931733069
            ,0.931982969
            ,0.932231526
            ,0.932477215
            ,0.93271802
            ,0.932954735
            ,0.933195052
            ,0.933433599
            ,0.93367068
            ,0.933903916
            ,0.934133368
            ,0.934361782
            ,0.934592576
            ,0.934817022
            ,0.935038661
            ,0.935261642
            ,0.935488408
            ,0.935720118
            ,0.935947494
            ,0.936164554
            ,0.936383324
            ,0.936603558
            ,0.9368172
            ,0.937036336
            ,0.937252297
            ,0.937458859
            ,0.937673477
            ,0.937883396
            ,0.938088981
            ,0.938296274
            ,0.938502592
            ,0.938700363
            ,0.938905642
            ,0.939110495
            ,0.939316873
            ,0.93952026
            ,0.939719863
            ,0.939921114
            ,0.940122364
            ,0.940320136
            ,0.940516992
            ,0.940713054
            ,0.940907407
            ,0.941097671
            ,0.941291536
            ,0.941481616
            ,0.941677251
            ,0.941865256
            ,0.942047767
            ,0.942233697
            ,0.94241981
            ,0.942601528
            ,0.942784832
            ,0.942969969
            ,0.943147841
            ,0.943324981
            ,0.943506759
            ,0.943686341
            ,0.943856034
            ,0.944033295
            ,0.944204575
            ,0.944372498
            ,0.944550981
            ,0.944723603
            ,0.944894944
            ,0.945064332
            ,0.945237687
            ,0.945407075
            ,0.945576645
            ,0.945742065
            ,0.945904006
            ,0.946065458
            ,0.94623564
            ,0.9463952
            ,0.946557934
            ,0.946720424
            ,0.946879435
            ,0.947034356
            ,0.947190742
            ,0.947352865
            ,0.94751389
            ,0.947673084
            ,0.94782886
            ,0.947986467
            ,0.948136993
            ,0.948289838
            ,0.948449704
            ,0.948599253
            ,0.948749779
            ,0.948898047
            ,0.949044544
            ,0.949197756
            ,0.949343643
            ,0.9494913
            ,0.949637981
            ,0.949782952
            ,0.949932684
            ,0.950080403
            ,0.950223054
            ,0.950368819
            ,0.950509823
            ,0.950652475
            ,0.950786947
            ,0.950926913
            ,0.951065109
            ,0.951206846
            ,0.951352855
            ,0.951489891
            ,0.951621372
            ,0.9517556
            ,0.951883236
            ,0.952021371
            ,0.952154623
            ,0.952290621
            ,0.952420088
            ,0.952552668
            ,0.95268653
            ,0.952820331
            ,0.952951446
            ,0.953079082
            ,0.953209953
            ,0.953333805
            ,0.953464798
            ,0.953588527
            ,0.953718726
            ,0.953847888
            ,0.953973327
            ,0.954095713
            ,0.954222128
            ,0.954342866
            ,0.954464215
            ,0.95458886
            ,0.954711185
            ,0.954829299
            ,0.954947534
            ,0.955069677
            ,0.955183822
            ,0.955298396
            ,0.955420355
            ,0.955537064
            ,0.955659695
            ,0.955779518
            ,0.955896593
            ,0.956014341
            ,0.956127327
            ,0.956245684
            ,0.956354947
            ,0.956467323
            ,0.956580309
            ,0.956692318
            ,0.95680561
            ,0.956921221
            ,0.957034573
            ,0.957145667
            ,0.957254014
            ,0.957365779
            ,0.957473333
            ,0.957574721
            ,0.957681115
            ,0.957787203
            ,1

        };

        static int generateTumorGermlineRNAReadDepth(Random random, int minReadDepth)
        {
            int result;

            do
            {
                var value = random.NextDouble();

                for (result = 0; result < tumorGermlineRNAReadDepthDistribution.Count() - 1; result++)
                {
                    if (tumorGermlineRNAReadDepthDistribution[result] > value)
                    {
                        break;
                    }
                }
            } while (result < minReadDepth);

            return result;
        }

        static double generateASE(Random random)
        {
            var randomValue = random.NextDouble();

            for (int i = 0; i < underlyingASEDistribution.Count(); i++)
            {
                if (underlyingASEDistribution[i] > randomValue)
                {
                    return (double)i / (underlyingASEDistribution.Count() - 1);
                }
            }
            return 1;
        } // generateASE

        static double generateMeasuredASE(Random random, double trueASE, int readCount)
        {
            int a = 0, b = 0;

            for (int i = 0; i < readCount; i++)
            {
                if (random.NextDouble() <= trueASE / 2 + .5)
                {
                    a++;
                } else
                {
                    b++;
                }
            }

            return (double)Math.Abs(a - b) / (a + b);
        }


        static int nIterations = 100000;
        static int maxMin = 250;
        static void Main(string[] args)
        {
            var outputFile = ASETools.CreateStreamWriterWithRetry(@"f:\temp\SimulatedASEDistribution.txt");
            var random = new System.Random();

            for (int readCount = 5; readCount <= maxMin; readCount += 5)
            {
                var differenceHistogram = new ASETools.Histogram("Difference -- Read Count " + readCount);
                var errorHistogram = new ASETools.Histogram("Error -- Read Count " + readCount);

                for (int i = 0; i < nIterations; i++)
                {
                    double realASE = generateASE(random);
                    var measuredASE1 = generateMeasuredASE(random, realASE, readCount);
                    var measuredASE2 = generateMeasuredASE(random, realASE, readCount);

                    errorHistogram.addValue(measuredASE1 - realASE);
                    errorHistogram.addValue(measuredASE2 - realASE);

                    differenceHistogram.addValue(Math.Abs(measuredASE1 - measuredASE2));
                }

                outputFile.WriteLine("Difference -- Read count " + readCount);
                outputFile.WriteLine(ASETools.HistogramResultLine.Header());
                differenceHistogram.ComputeHistogram(0, 1, 0.01).ToList().ForEach(x => outputFile.WriteLine(x.ToString()));

                outputFile.WriteLine("Error -- Read Count " + readCount);
                outputFile.WriteLine(ASETools.HistogramResultLine.Header());
                errorHistogram.ComputeHistogram(-1, 1, 0.01).ToList().ForEach(x => outputFile.WriteLine(x.ToString()));
                Console.Write(".");
            }

            var realASEDistribution = new Dictionary<int, Dictionary<int, ASETools.RunningMeanAndStdDev>>();   // Maps read depth -> measured ASE times 100 -> real ASEs.  Only done for min read count 10.

            for (int minReadCount = 10; minReadCount <= maxMin; minReadCount += 5)
            {
                var differenceHistogram = new ASETools.Histogram();
                var errorHistogram = new ASETools.Histogram();

                for (int i = 0; i < nIterations; i++)
                {
                    var readCount1 = generateTumorGermlineRNAReadDepth(random, minReadCount);
                    var readCount2 = generateTumorGermlineRNAReadDepth(random, minReadCount);
                    double realASE = generateASE(random);

                    var measuredASE1 = aseCorrection.getCorrectedASE( generateMeasuredASE(random, realASE, readCount1), readCount1);
                    var measuredASE2 = aseCorrection.getCorrectedASE(generateMeasuredASE(random, realASE, readCount2), readCount2);

                    errorHistogram.addValue(measuredASE1 - realASE);
                    errorHistogram.addValue(measuredASE2 - realASE);

                    differenceHistogram.addValue(Math.Abs(measuredASE1 - measuredASE2));

                    if (minReadCount == 10) // We only compute the back matrix for the smallest min read depth
                    {
                        if (!realASEDistribution.ContainsKey(readCount1))
                        {
                            realASEDistribution.Add(readCount1, new Dictionary<int, ASETools.RunningMeanAndStdDev>());
                        }

                        int roundedMeasuredASETimesOneHundred = (int)Math.Round(measuredASE1 * 100);

                        if (!realASEDistribution[readCount1].ContainsKey(roundedMeasuredASETimesOneHundred))
                        {
                            realASEDistribution[readCount1].Add(roundedMeasuredASETimesOneHundred, new ASETools.RunningMeanAndStdDev());
                        }

                        realASEDistribution[readCount1][roundedMeasuredASETimesOneHundred].addValue(realASE);
                    }
                }

                outputFile.WriteLine("Difference -- Read count selected from real distribution with min " + minReadCount);
                outputFile.WriteLine(ASETools.HistogramResultLine.Header());
                differenceHistogram.ComputeHistogram(0, 1, 0.01).ToList().ForEach(x => outputFile.WriteLine(x.ToString()));

                outputFile.WriteLine("Error -- Read count selected from real distribution with min " + minReadCount);
                outputFile.WriteLine(ASETools.HistogramResultLine.Header());
                errorHistogram.ComputeHistogram(-1, 1, 0.01).ToList().ForEach(x => outputFile.WriteLine(x.ToString()));

                Console.Write(".");
            }

            outputFile.WriteLine("Distribution of real ASE given read depth and measured ASE");
            outputFile.Write("Read Depth");
            for (int ase = 0; ase <= 100; ase++)
            {
                outputFile.Write("\t" + (double)ase/100.0);
            }
            outputFile.WriteLine();
            for (int readDepth = 10; readDepth <= 250; readDepth++)
            {
                outputFile.Write(readDepth);
                for (int ase = 0; ase <= 100; ase++)
                {
                    if (realASEDistribution.ContainsKey(readDepth) && realASEDistribution[readDepth].ContainsKey(ase))
                    {
                        outputFile.Write("\t" + realASEDistribution[readDepth][ase].getMeanAndStdDev().mean);
                    } else
                    {
                        outputFile.Write("\t*");
                    }
                }
                outputFile.WriteLine();
            }

            outputFile.WriteLine("Distribution of real ASE minus measured ASE given read depth and measured ASE");
            outputFile.Write("Read Depth");
            for (int ase = 0; ase <= 100; ase++)
            {
                outputFile.Write("\t" + (double)ase / 100.0);
            }
            outputFile.WriteLine();
            for (int readDepth = 10; readDepth <= 250; readDepth++)
            {
                outputFile.Write(readDepth);
                for (int ase = 0; ase <= 100; ase++)
                {
                    if (realASEDistribution.ContainsKey(readDepth) && realASEDistribution[readDepth].ContainsKey(ase))
                    {
                        outputFile.Write("\t" + (realASEDistribution[readDepth][ase].getMeanAndStdDev().mean - (double)ase / 100.0));
                    }
                    else
                    {
                        outputFile.Write("\t*");
                    }
                }
                outputFile.WriteLine();
            }



            outputFile.WriteLine("**done**");
            outputFile.Close();

        }
    }
}
