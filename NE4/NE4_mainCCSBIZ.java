package paper1.NE0602.NE4;

import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class NE4_mainCCSBIZ {
	public static void main(String[] args) {
        long startTime = System.currentTimeMillis();
        int repeat = 1;
        int correct = 0;
        double sampleSize = 0d;
        long seed= 100;
        Random R = new Random(seed);

        ArrayList<int[]> listALT = new ArrayList<>();
        int altcount = 0;
        int budget1 = 22;
        int budget2 = 13;
        int ac = 1;
        int bc = 1;
        int cc = 2;
        int dc = 2;
        int ec = 1;
        int fc = 1;

        double[][] bound = new double[][]{
                {1.0, 10}, //room
                {1.0, 10}, //cubicle
                {1.0, 6}, //resuscitation
                {1.0, (double)budget2/dc}, //labtest
                {1.0, (double)budget2/ec}, //registration
                {1.0, (double)budget2/fc}, //triage
        };
        for (int a = 1; a <= bound[0][1]; a++){
            for (int b = 1; b <= bound[1][1]; b++){
                for (int c = 1; c <= bound[2][1]; c++){
                    for (int d = 1; d <= bound[3][1]; d++){
                        for (int e = 1; e <= bound[4][1]; e++){
                            for (int f = 1; f <= bound[5][1]; f++){
                                if (a*ac+b*bc+c*cc==budget1 && d*dc+e*ec+f*fc==budget2){
                                    int[] alt = new int[6];
                                    alt[0] = a;
                                    alt[1] = b;
                                    alt[2] = c;
                                    alt[3] = d;
                                    alt[4] = e;
                                    alt[5] = f;
                                    listALT.add(alt);
                                    altcount+=1;
                                }
                            }
                        }
                    }
                }
            }
        }


        int k = altcount;
        int d = 3;

        double feasibleThreshold = 0.5;
        double alpha = 0.05;

        ArrayList<Double> gammas = new ArrayList<Double>();
        ArrayList<Double> ssinterval= new ArrayList<>();


        for (int i = 0; i < d; i++) {
            gammas.add(0.1);
        }
        double izpp = 2;
        ArrayList<Double> izsp = new ArrayList<Double>();

        izsp.add(0.03);
        izsp.add(0.03);
        izsp.add(0.03);
        mnCalForCCSBIZ n0cal = new mnCalForCCSBIZ(alpha, d, k, gammas, izsp);
        n0cal.run();
        int n00 = n0cal.getn0();
        ArrayList<Integer> mbetan0 = n0cal.getmbetan0();
        double alpha2 = n0cal.getAlpha2();

        System.out.println("n00:"+n00);
        System.out.println(mbetan0);


        int rows = repeat;
        int cols = 6;
        double[][] data = new double[rows][cols];


        for (int count = 0; count < repeat; count++){
            System.out.println(count);
            long tempSeed = R.nextLong();
            NE4_CCSBIZ y1 = new NE4_CCSBIZ(alpha, d, gammas, feasibleThreshold, listALT, tempSeed, izpp, izsp, n00, mbetan0, alpha2);
            y1.run();
//            if(y1.getBestID()==k/2-1) {
//                correct +=1;
//            }
            data[count][0] = count;
            data[count][1] = y1.getBestID();
            data[count][2] = y1.getTotalSampleSize();
            data[count][3] = seed;
            data[count][4] = izpp;
            data[count][5] = izsp.get(0);

            if(y1.getBestID()==680) {
                correct +=1;
            }
            sampleSize += y1.getTotalSampleSize();
            ssinterval.add(y1.getTotalSampleSize());
        }
        double std = sd(ssinterval, sampleSize/repeat);
        int halfinterval = (int)(1.96*std/Math.sqrt(repeat));


        String filePath = "C:/Users/BRZ/Desktop/NE/NE4_0517/NE5_0525/NE4Hong"+seed+".csv";  // 结果输出路径
        try (FileWriter writer = new FileWriter(filePath)) {
            StringBuilder header = new StringBuilder();

            header.append("index");
            header.append(",");
            header.append("bestalt");
            header.append(",");
            header.append("samplesize");
            header.append(",");
            header.append("seed");
            header.append(",");
            header.append("izpp");
            header.append(",");
            header.append("izsp");


            writer.write(header.toString());
            writer.write("\n");

            for (int i = 0; i < rows; i++) {
                StringBuilder row = new StringBuilder();
                for (int j = 0; j < cols; j++) {
                    row.append(data[i][j]);
                    if (j < cols - 1) {
                        row.append(",");
                    }
                }
                writer.write(row.toString());
                writer.write("\n");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }


        long endTime = System.currentTimeMillis();
        long runTime = endTime - startTime;
	}

    public static double sd(ArrayList<Double> al, double average) {
        double sum = 0;
        for (int i = 0; i < al.size(); i++) {
            sum = sum + (al.get(i) - average) * (al.get(i) - average);
        }
        double standardDeviation = Math.sqrt(sum/(al.size()-1));
        return standardDeviation;
    }

    public static double inverseCumulativeProbabilityforSTD(double p){
        NormalDistribution normalDistribution = new NormalDistribution();
        return normalDistribution.inverseCumulativeProbability(p);
    }
}
