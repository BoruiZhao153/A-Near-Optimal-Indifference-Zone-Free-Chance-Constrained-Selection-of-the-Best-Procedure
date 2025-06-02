package paper1.NE0602.NE2;

import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.ArrayList;
import java.util.Random;

public class NE2_mainCCSBIZ {
	public static void main(String[] args) {
        long startTime = System.currentTimeMillis();
        int repeat = 1000;
        int correct = 0;
        double sampleSize = 0d;
        long seed= 1;
        Random R = new Random(seed);

        int k = 501;
        int d = 6;
        int b = (int)Math.floor(1.0*(k+1)/2); //251
//        System.out.println(b+" "+(b+(k-b)/3)+" "+(b+2*(k-b)/3));
//        System.out.println((double)(k-b)/3);
        int n0 = 10;
        double feasibleThreshold = 0.5;
        double deltax = 1/Math.sqrt(n0);
        double deltay = 0.02;
        double alpha = 0.05;
//        double theta = 0.5;
//        double cThreshold=0;
        ArrayList<Double> gammas = new ArrayList<Double>();
        ArrayList<Double> ssinterval= new ArrayList<>();

        int caseIndex = 3;


        for (int i = 0; i < d; i++) {
            gammas.add(0.1);
        }

        double izpp = 2*deltax;
        ArrayList<Double> izsp = new ArrayList<Double>();
        for (int i = 0; i < d; i++) {
            izsp.add(2*deltay);
        }

        mnCalForHong_2 n0cal = new mnCalForHong_2(alpha, d, k, gammas, izsp);
        n0cal.run();
        int n00 = n0cal.getn0();
        ArrayList<Integer> mbetan0 = n0cal.getmbetan0();
        double alpha2 = n0cal.getAlpha2();
        System.out.println("n00:"+n00);
        System.out.println(mbetan0);
        ArrayList<Double> gssinterval= new ArrayList<>();
        for (int count = 0; count < repeat; count++){
            System.out.println(count);

            ArrayList<Double> ppmu = new ArrayList<Double>();
            ArrayList<Double> ppsigma2 = new ArrayList<Double>();
            ArrayList<ArrayList<Double>> spberp = new ArrayList<ArrayList<Double>>();
            for (int i = 0; i < k; i++) {
                spberp.add(new ArrayList<>());
            }

            ArrayList<int[]> violateConstrainIndex = new ArrayList<int[]>();
            for(int i = 0; i < k; i++) {
                if(i == 0){
                    violateConstrainIndex.add(new int[]{0, 0, 0, 0, 0, 0});
                }else if(i < b){
                    violateConstrainIndex.add(new int[]{0, 0, 0, 0, 0, 0});
                }else{
                    if (caseIndex==1){
                        violateConstrainIndex.add(new int[]{1, 1, 0, 0, 0, 0});
                    }else if(caseIndex==2){
                        violateConstrainIndex.add(new int[]{1, 1, 1, 1, 0, 0});
                    }else if(caseIndex==3){
                        violateConstrainIndex.add(new int[]{1, 1, 1, 1, 1, 1});
                    }
                }
//                System.out.println(i+" "+ Arrays.toString(violateConstrainIndex.get(i)));
            }

            for(int i = 0; i < k; i++) {
                if(i == 0){
                    ppmu.add(3*deltax);
                }else if(i < b){
                    ppmu.add(R.nextDouble()*2*deltax);
                }else{
                    ppmu.add(R.nextDouble()*5*deltax);
                }
                ppsigma2.add(100.0); //////////////////////

                for (int j = 0; j < d; j++){
                    if (violateConstrainIndex.get(i)[j]==1){
                        spberp.get(i).add(1-gammas.get(j)-R.nextDouble()*0.4);
                    }else{
                        if (i==0){
                            spberp.get(i).add(1-gammas.get(j)+deltay);
                        }else{
                            spberp.get(i).add(1-gammas.get(j)+R.nextDouble()*0.1);
                        }
                    }
                }
            }


            long tempSeed = R.nextLong();
            NE2_CCSBIZ y1 = new NE2_CCSBIZ(alpha, d, gammas, feasibleThreshold, ppmu, ppsigma2, spberp, tempSeed, izpp, izsp, n00, mbetan0, alpha2);
            y1.run();
//            if(y1.getBestID()==k/2-1) {
//                correct +=1;
//            }
            if(y1.getBestID()==0) {
                correct +=1;
            }

            sampleSize += y1.getTotalSampleSize();
            ssinterval.add(y1.getTotalSampleSize());

        }
        double std = sd(ssinterval, sampleSize/repeat);
        int halfinterval = (int)(1.96*std/Math.sqrt(repeat));
        System.out.println("ba样本量0.95的置信区间为:  " + sampleSize/repeat + "±" + halfinterval);

        System.out.println(correct*1d/repeat+" "+sampleSize/repeat);

        long endTime = System.currentTimeMillis();
        long runTime = endTime - startTime;
        System.out.println(izpp+" "+izsp+" caseindex " +caseIndex);
        System.out.println("代码的运行时间为： "+runTime+ " 毫秒");
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
