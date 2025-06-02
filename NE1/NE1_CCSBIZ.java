package paper1.NE0602.NE1;

import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.ArrayList;
import java.util.Random;

public class NE1_CCSBIZ {
    private int k;
    private double alpha;
    private double izpp;  //indifference-zone parameter
    private double izsp;  // feasibility tolerance parameter
    private double y;
    private double feasibleThreshold;
    private int bestID = -1;
    private double totalSampleSize=0;
    private ArrayList<Double> ppmu = new ArrayList<Double>();  //mu of the primary performance(standard)
    private ArrayList<Double> ppsigma2 = new ArrayList<Double>();  //sigma2 of primary performance(standard)
    private ArrayList<Double> spberp = new ArrayList<Double>();
    private ArrayList<int[]> listALT;
    private Random R= new Random();

    public NE1_CCSBIZ(double alpha, double izpp, double izsp, double y, double feasibleThreshold, ArrayList<Double> ppmu
            , ArrayList<Double> ppsigma2, ArrayList<Double> spberp, long seed){
        this.alpha=alpha;
        this.izpp=izpp;
        this.izsp=izsp;
        this.y=y;
        this.feasibleThreshold=feasibleThreshold;
        this.ppmu.addAll(ppmu);
        this.ppsigma2.addAll(ppsigma2);
        this.spberp.addAll(spberp);
        R.setSeed(seed);
    }
    public NE1_CCSBIZ(){
    }

    public int getBestID() {
        return bestID;
    }

    public double getTotalSampleSize() {
        return totalSampleSize;
    }

    public void run() {
        k  = ppmu.size();
        double alpha1 = alpha/k;
        double alpha2 = alpha/k;
        long startTime1 = System.currentTimeMillis();
        ArrayList<Integer> I = new ArrayList<Integer>();
        for(int i = 0 ; i < k ; i++) {
            I.add(i);
        }

        double n0e = inverseCumulativeProbabilityforSTD(1-alpha1)*inverseCumulativeProbabilityforSTD(1-alpha1)/(izsp*izsp)
                *(Math.sqrt((y-izsp)*(1-y+izsp))+Math.sqrt(y*(1-y)))*(Math.sqrt((y-izsp)*(1-y+izsp))+Math.sqrt(y*(1-y)));
        int n0 = (int)Math.ceil(n0e);
        System.out.println("n0e: "+n0e+" n0: "+n0);
        double malpha1n0e = n0*y - inverseCumulativeProbabilityforSTD(1-alpha1)*Math.sqrt(n0*y*(1-y));
        int malpha1n0 = (int)Math.floor(malpha1n0e);
        System.out.println("malpha1n0e: "+malpha1n0e+" malpha1n0: "+malpha1n0);
        double h2 = (n0 - 1)*(Math.pow(2*alpha2,-2.0/(n0-1))-1);

        // Feasibility check
        int fcellnum = 0;
        int ocellnum = 0;
        ArrayList<Integer> F = new ArrayList<Integer>(); //I feasible
        ArrayList<Integer> Fold = new ArrayList<Integer>(); //I feasible
        String[] altFeasible = new String[k];
        for (int i = 0; i < altFeasible.length; i++) {
            altFeasible[i] = "unknown";
        }
        int i1 = 0;
        double numsamplefeasibleell = 0.0;
        double numsamplefeasiblepass = 0.0;
        double numsampleoptimalonly = 0.0;
        double numaltinoptimalwaste = 0.0;
        double[][] sampleinfea =  new double[k][n0];
        while (i1<k){
            int tau = 0;
            int Ztau = 0;
            while (altFeasible[i1].equals("unknown")){
                tau = tau + 1;

                sampleinfea[i1][tau-1] = R.nextGaussian() * Math.sqrt(ppsigma2.get(I.get(i1)))+ppmu.get(I.get(i1));
                double samplesp = bernoullisample(spberp.get(I.get(i1)));

                if (samplesp < 0.5){
                    Ztau = Ztau + 1;
                }
                if (i1 == 0){
//                    System.out.println(samplesp1);
//                    System.out.println((double)Ztau/tau);
                }
                if(Ztau >= malpha1n0 + 1){
//                    System.out.println("Ztau: "+Ztau+" tau: "+tau);
                    altFeasible[i1] = "infeasible";
                    fcellnum = fcellnum + 1;
//                    System.out.println(tau+"时长通过可行性分析淘汰了"+I.get(i1)+" 其Ztau均值为："+(double)Ztau/tau+" 可行性分析已淘汰" + fcellnum);
                    totalSampleSize = totalSampleSize + tau;
                    numsamplefeasibleell = numsamplefeasibleell + tau;
                }else if(tau == n0){
                    ////////////////////////////////////////////////////////////////////////
                    altFeasible[i1] = "feasible";
                    numsamplefeasiblepass = numsamplefeasiblepass+tau;
//                    System.out.println(tau+"时长"+I.get(i1)+"通过了可行性分析， 其Ztau均值为："+(double)Ztau/tau);
                    F.add(i1);
                }
            }
            i1 = i1 + 1;
//            System.out.println(i1);
//            System.out.println(Arrays.deepToString(altFeasible));
//            System.out.println(F);
        }
        System.out.println("可行方案集F的大小为： "+F.size()+" 可行性分析进行淘汰花费的样本数为："+numsamplefeasibleell+" 通过可行性分析花费的样本数为："+numsamplefeasiblepass);
        long endTime2 = System.currentTimeMillis();
        long runTime2 = endTime2 - startTime1;
        System.out.println("可行性分析后方案集F的大小为： "+F.size()+" 样本量为："+totalSampleSize+" 花费时间为："+runTime2);

        //Selecting the Best
        double[] tempX = new double[F.size()];
        double[][] Xij =  new double[k][k];
        double[][] Xij2 = new double[k][k];
        double[][] Sij2 = new double[k][k];
        for(int ell = 0; ell < n0; ell++) {
            for(int i = 0; i < F.size(); i++) {
                tempX[i]= sampleinfea[F.get(i)][ell];
            }

            for(int i = 0; i < F.size(); i++) {
                Xij[F.get(i)][F.get(i)] = Xij[F.get(i)][F.get(i)] + tempX[i];
                Xij2[F.get(i)][F.get(i)] = Xij2[F.get(i)][F.get(i)] + tempX[i] * tempX[i];

                for(int j = i + 1; j < F.size(); j++) {
                    Xij[F.get(i)][F.get(j)] = Xij[F.get(i)][F.get(j)]+tempX[i]-tempX[j];
                    Xij[F.get(j)][F.get(i)] = Xij[F.get(j)][F.get(i)]+tempX[j]-tempX[i];
                    Xij2[F.get(i)][F.get(j)] = Xij2[F.get(i)][F.get(j)] +(tempX[i]-tempX[j])*(tempX[i]-tempX[j]);
                    Xij2[F.get(j)][F.get(i)] = Xij2[F.get(i)][F.get(j)];
                }
            }
        }

        for(int i = 0; i < F.size(); i++) {
            for(int j = i + 1; j < F.size(); j++) {
                Sij2[F.get(i)][F.get(j)] = (Xij2[F.get(i)][F.get(j)] - Xij[F.get(i)][F.get(j)]* Xij[F.get(i)][F.get(j)]/n0)/(n0-1);
                Sij2[F.get(j)][F.get(i)] = Sij2[F.get(i)][F.get(j)];
            }
        }

        int r = n0;

        while(F.size()>1) {
            for(int i  = 0;  i < F.size(); i++) {
                for (int j = 0; j < F.size(); j++) {
                    if(i!=j) {
                        double Zij = Xij[F.get(i)][F.get(j)];
                        double gt = h2*Sij2[F.get(i)][F.get(j)]/(2.0*izpp)-izpp*r/2d;
                        if (gt < 0){
//                            System.out.println("gt<0000000000000000000000000000000000000000000000000000000000000000000000");
                            gt = 0;
                        }
                        if(Zij > gt) {
                            totalSampleSize = totalSampleSize + r;
                            numsampleoptimalonly = numsampleoptimalonly + r - n0;
                            ocellnum = ocellnum + 1;
                            if (r == n0){
                                numaltinoptimalwaste = numaltinoptimalwaste + 1;
                            }
//                            System.out.println(r+"时"+F.get(i)+"凭借pp为"+Xij[F.get(i)][F.get(i)]/r+"的表现"+"通过最优性分析去掉了"+F.get(j)+" 已经通过可行性分析淘汰"+fcellnum
//                                    +" 已通过最优性分析淘汰"+ocellnum+" 总共淘汰了"+(fcellnum+ocellnum)+" 其pp为："+Xij[F.get(j)][F.get(j)]/r);

                            F.remove(j);
                            if(j < i) {
                                i--;
                            }
                            j--;
                        }
                    }
                }
            }
            if(F.size()==1) {
                totalSampleSize = totalSampleSize + r;
                numsampleoptimalonly = numsampleoptimalonly + r - n0;
            }
            r = r+1;
            double[] tempX2 = new double[F.size()];
            for(int i = 0; i < F.size(); i++) {
                tempX2[i]= R.nextGaussian() * Math.sqrt(ppsigma2.get(F.get(i)))+ppmu.get(F.get(i));
            }

            for(int i = 0; i < F.size(); i++) {
                Xij[F.get(i)][F.get(i)] = Xij[F.get(i)][F.get(i)] + tempX2[i];
            }

            for(int i = 0; i < F.size(); i++) {
                for(int j = i + 1; j < F.size(); j++) {
                    Xij[F.get(i)][F.get(j)] = Xij[F.get(i)][F.get(j)]+tempX2[i]-tempX2[j];
                    Xij[F.get(j)][F.get(i)] = Xij[F.get(j)][F.get(i)]+tempX2[j]-tempX2[i];
                }
            }
        }
        System.out.println("可行性分析进行淘汰花费的样本数为："+numsamplefeasibleell+" 通过可行性分析花费的样本数为："+numsamplefeasiblepass
                +" 最优性分析使用的额外样本数为："+numsampleoptimalonly+" n0对于最优性分析已经过剩的方案数为："+numaltinoptimalwaste+"占此阶段总淘汰方案数的比例为："+(double)numaltinoptimalwaste/ocellnum);
        bestID = F.get(0);
        long endTime1 = System.currentTimeMillis();
        long runTime1 = endTime1 - startTime1;
        System.out.println("最优方案为："+bestID+" 样本量为："+totalSampleSize+" 花费时间为："+runTime1);


    }

    public double inverseCumulativeProbabilityforSTD(double p){
        NormalDistribution normalDistribution = new NormalDistribution();
        return normalDistribution.inverseCumulativeProbability(p);
    }

    public double bernoullisample(double p){
//        System.out.println(R.nextDouble());
        boolean result = R.nextDouble() < p;

        if (result) {
            return 1;
        } else {
            return 0;
        }
    }
}
