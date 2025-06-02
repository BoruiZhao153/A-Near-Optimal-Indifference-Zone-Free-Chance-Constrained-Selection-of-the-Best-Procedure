package paper1.NE0602.NE2;

import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.ArrayList;
import java.util.Random;

public class NE2_mainP2 {
	public static void main(String[] args) {
        long startTime = System.currentTimeMillis();
        int repeat = 1000;
        int correct = 0;
        double sampleSize = 0d;
        long seed= 1;
        Random R = new Random(seed);

        int k = 501;
        int d = 6;
        int b = (int)Math.floor(1.0*(k+1)/2);
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

        double izpp = 2*deltax;
        ArrayList<Double> izsp = new ArrayList<Double>();
        for (int i = 0; i < d; i++) {
            izsp.add(2*deltay);
        }
//        CCSP1 y1 = new CCSP1();
        int goodalternativecorrect = 0;
        int goodalternativefalse = 0;
        int findgoodalternativecount = 0;
        double goodalternativesamplesize = 0.0;
        ArrayList<Double> gssinterval= new ArrayList<>();
        System.out.println(caseIndex);
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
            NE2_P2 y1 = new NE2_P2(alpha, d, gammas, ppmu, ppsigma2, spberp, tempSeed, izpp, izsp);
            y1.run();
//            if(y1.getBestID()==k/2-1) {
//                correct +=1;
//            }
            if(y1.getBestID()==0) {
                correct +=1;
            }
            if(y1.getGoodIDCorrectOrNot()==1){
                goodalternativecorrect+=1;
                findgoodalternativecount+=1;
                goodalternativesamplesize+=y1.getGoodAlternativeSampleSize();
                gssinterval.add(y1.getGoodAlternativeSampleSize());
            }else if(y1.getGoodIDCorrectOrNot()==0){
                goodalternativefalse+=1;
                findgoodalternativecount+=1;
                goodalternativesamplesize+=y1.getGoodAlternativeSampleSize();
                gssinterval.add(y1.getGoodAlternativeSampleSize());
            }
            sampleSize += y1.getTotalSampleSize();
            ssinterval.add(y1.getTotalSampleSize());
        }
        double std = sd(ssinterval, sampleSize/repeat);
        int halfinterval = (int)(1.96*std/Math.sqrt(repeat));
        System.out.println("ba样本量0.95的置信区间为:  " + sampleSize/repeat + "±" + halfinterval);

        System.out.println(correct*1d/repeat+" "+sampleSize/repeat);

        double gstd = sd(gssinterval, goodalternativesamplesize/findgoodalternativecount);
        int ghalfinterval = (int)(1.96*gstd/Math.sqrt(findgoodalternativecount));
        System.out.println("ga样本量0.95的置信区间为:  " + goodalternativesamplesize/findgoodalternativecount + "±" + ghalfinterval);

        System.out.println("选出ga且选择正确的概率为"+goodalternativecorrect*1d/findgoodalternativecount+" 选择出正确goodalternative的概率为"+goodalternativecorrect*1d/repeat+" 选择出错误goodalternative的概率为"+goodalternativefalse*1d/repeat
                +" 未选择出goodalternative的概率为"+(1-goodalternativecorrect*1d/repeat-goodalternativefalse*1d/repeat));

        long endTime = System.currentTimeMillis();
        long runTime = endTime - startTime;
        System.out.println("caseindex " +caseIndex);
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
