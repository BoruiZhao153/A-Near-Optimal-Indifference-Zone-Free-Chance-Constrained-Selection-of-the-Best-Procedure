package paper1.NE0602.NE3;

import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.ArrayList;
import java.util.Random;

public class NE3_mainP2 {
	public static void main(String[] args) {
        long startTime = System.currentTimeMillis();
        int repeat = 1000;
        int correct = 0;
        double sampleSize = 0d;
        long seed= 1;
        Random R = new Random(seed);

        int k = 504;
        int d = 6;
        int b = (int)Math.floor(1.0*(k+1)/2); //252
        int n0 = 10;

        double deltax = 1/Math.sqrt(n0);
        double deltay = 0.02;
        double alpha = 0.05;

        ArrayList<Double> gammas = new ArrayList<Double>();
        ArrayList<Double> ssinterval= new ArrayList<>();

        int caseIndex = 3;

        for (int i = 0; i < d; i++) {
            gammas.add(0.1);
        }

        double izpp = 1*deltax;
        ArrayList<Double> izsp = new ArrayList<Double>();
        for (int i = 0; i < d; i++) {
            izsp.add(1*deltay);
        }

        double goodalternativesamplesize = 0.0;
        ArrayList<Double> gssinterval= new ArrayList<>();

        double BG = 0; //BEST AS GOOD ALT
        double NBG = 0; //GOOD ALT BUT NOT BEST
        double NG = 0; //NOT GOOD

        for (int count = 0; count < repeat; count++){
            System.out.println("count: "+count);

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
            }
            for(int i = 0; i < k; i++) {
                if(i == 0){
                    ppmu.add(3*deltax);
                }else if(i == 1){
                    ppmu.add(2.9*deltax);
                }else if(i < b){
                    ppmu.add(R.nextDouble()*2*deltax);
                }else if(i == b){
                    ppmu.add(2.9*deltax);
                }else if(i == b + 1){
                    ppmu.add(3.1*deltax);
                }else{
                    ppmu.add(R.nextDouble()*5*deltax);
                }
                ppsigma2.add(100.0); //////////////////////

                for (int j = 0; j < d; j++){
                    if (violateConstrainIndex.get(i)[j]==1){
                        if (i==b || i==b+1){
                            spberp.get(i).add(1-gammas.get(j)-0.1*deltay);
                        }else {
                            spberp.get(i).add(1-gammas.get(j)-deltay-R.nextDouble()*deltay*2);
                        }
                    }else{
                        if (i==0){
                            spberp.get(i).add(1-gammas.get(j)+deltay);
                        }else if (i==1){
                            spberp.get(i).add(1-gammas.get(j)+0.1*deltay);
                        }else {
                            spberp.get(i).add(1-gammas.get(j)+deltay+R.nextDouble()*deltay*2);
                        }
                    }
                }
                if (spberp.get(i).size()>d){
                    System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                }
            }
            int numGoodAlt = 0;
            for (int i = 0; i < k; i++ ){
                boolean goodAltOrNot = true;
                for (int j = 0; j < d; j++){
                    if (spberp.get(i).get(j) < 1 - gammas.get(j) - izsp.get(j)){
                        goodAltOrNot = false;
                    }
                }
                if (ppmu.get(0) > ppmu.get(i)+izpp){
                    goodAltOrNot = false;
                }
                if (goodAltOrNot){
                    numGoodAlt+=1;
                }
            }
            System.out.println((double)izpp/deltax+" "+(double)izsp.getFirst()/deltay+" numgoodalt: "+numGoodAlt);

            long tempSeed = R.nextLong();
            NE3_P2 y1 = new NE3_P2(alpha, d, gammas, ppmu, ppsigma2, spberp, tempSeed, izpp, izsp);
            y1.run();

            if(y1.getBestID()==0) {
                correct +=1;
            }

            if (y1.getGoodID()==0){
                BG+=1;
                goodalternativesamplesize+=y1.getGoodAlternativeSampleSize();
                gssinterval.add(y1.getGoodAlternativeSampleSize());
            }else if (y1.getGoodID()!=0 && y1.getGoodIDCorrectOrNot()==1){
                NBG+=1;
                goodalternativesamplesize+=y1.getGoodAlternativeSampleSize();
                gssinterval.add(y1.getGoodAlternativeSampleSize());
            }else if(y1.getGoodIDCorrectOrNot()==0){
                NG+=1;
                goodalternativesamplesize+=y1.getGoodAlternativeSampleSize();
                gssinterval.add(y1.getGoodAlternativeSampleSize());
            }

            sampleSize += y1.getTotalSampleSize();
            ssinterval.add(y1.getTotalSampleSize());
        }
        double std = sd(ssinterval, sampleSize/repeat);
        int halfinterval = (int)(1.96*std/Math.sqrt(repeat));
//        System.out.println("ba样本量0.95的置信区间为:  " + sampleSize/repeat + " ± " + halfinterval);
//        System.out.println(correct*1d/repeat+" "+sampleSize/repeat+" ± " + halfinterval);

        double gstd = sd(gssinterval, goodalternativesamplesize/repeat);
        int ghalfinterval = (int)(1.96*gstd/Math.sqrt(repeat));
//        System.out.println("ga样本量0.95的置信区间为:  " + goodalternativesamplesize/repeat + " ± " + ghalfinterval);
//
//
//        System.out.println("BG "+BG+" "+BG/repeat);
//        System.out.println("NBG "+NBG+" "+NBG/repeat);
//        System.out.println("BG "+NG+" "+NG/repeat);

        long endTime = System.currentTimeMillis();
        long runTime = endTime - startTime;
        System.out.println(izpp+" "+izsp);
//        System.out.println("代码的运行时间为： "+runTime+ " 毫秒");
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
