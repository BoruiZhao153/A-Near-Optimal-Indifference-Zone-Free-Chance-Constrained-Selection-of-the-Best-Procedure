package paper1.NE0602.NE1;

import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.ArrayList;
import java.util.Random;

public class NE1_mainP1 {
	public static void main(String[] args) {
        long startTime = System.currentTimeMillis();
        int repeat = 1000;
        int correct = 0;
        double sampleSize = 0d;
        long seed= 1;
        Random R = new Random(seed);

        int k = 501;
        int b = (int)Math.floor(1.0*(k+1)/2); //251
        int n0 = 10;
        double feasibleThreshold = 0;
        double gamma = 0.1; //violation error
        double deltax = 1/Math.sqrt(n0);
        double deltay = 0.02;

        double alpha = 0.05;

        ArrayList<Double> ssinterval= new ArrayList<>();


//        CCSP1 y1 = new CCSP1();
        System.out.println("k: "+k+" b-1 "+(b-1));
        for (int count = 0; count < repeat; count++){
            System.out.println("di"+count);

            ArrayList<Double> ppmu = new ArrayList<Double>();
            ArrayList<Double> ppsigma2 = new ArrayList<Double>();
            ArrayList<Double> spberp = new ArrayList<Double>();

            for(int i = 0; i < k; i++) {
                if(i == 0){
                    ppmu.add(3*deltax);
                    spberp.add(1-gamma+deltay);
                }else if(i < b){
                    ppmu.add(R.nextDouble()*2*deltax);
                    spberp.add(1-gamma+0.1/(b-1)*(i));
//                    System.out.println(1-gamma+0.1/250*(i));
                }else {
                    ppmu.add(R.nextDouble()*5*deltax);
                    spberp.add(1-gamma-0.4/(b-1)*(i-b+1));
//                    System.out.println(1-gamma-0.4/250*(i-b+1));
                }
                ppsigma2.add(100.0);
            }

            long tempSeed = R.nextLong();
            NE1_P1 y1 = new NE1_P1(alpha, gamma, feasibleThreshold, ppmu, ppsigma2, spberp, tempSeed);
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
        System.out.println("样本量0.95的置信区间为:  " + sampleSize/repeat + "±" + halfinterval);

        System.out.println(correct*1d/repeat+" "+sampleSize/repeat);
        long endTime = System.currentTimeMillis();
        long runTime = endTime - startTime;
        System.out.println(k);
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
