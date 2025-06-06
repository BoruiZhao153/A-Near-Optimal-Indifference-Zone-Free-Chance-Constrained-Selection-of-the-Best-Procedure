package paper1.NE0602.NE4;

import java.util.*;

public class NE4_ED {
//    int seed;
    int[] alt = new int[6];
    private double averageWTForC4C5;

    private double satisfyC1OrNot = 0;
    private double satisfyC2OrNot = 0;
    private double satisfyC3OrNot = 0;
    private long sortingTime = 0;
    static Random R = new Random();

    public NE4_ED(long seed){
        R.setSeed(seed);
    }
//    public EmergencyDepartmentSimulation1(int seed){
//        R.setSeed((long)seed);
//    }
    public void setAlt(int[] alt) {
        this.alt[0]=alt[0];this.alt[1]=alt[1];this.alt[2]=alt[2];this.alt[3]=alt[3];
        this.alt[4]=alt[4];this.alt[5]=alt[5];
    }

    public double getSatisfyC1OrNot(){return satisfyC1OrNot;}
    public double getSatisfyC2OrNot(){return satisfyC2OrNot;}
    public double getSatisfyC3OrNot(){return satisfyC3OrNot;}

    public double getAverageWTForC4C5(){return averageWTForC4C5;}
    public double getSortingTime(){return sortingTime;}



    public void run() {

        long startTime = System.currentTimeMillis();
        sortingTime = 0;
        ArrayList<Integer> patientId = new ArrayList<>();
        ArrayList<String> patientType = new ArrayList<>();
        ArrayList<String> currentNode = new ArrayList<>();
        ArrayList<String> nextNode = new ArrayList<>();
        ArrayList<Double> nextNodeTime = new ArrayList<>();
        ArrayList<Double> arrivalTimeCurrentNode = new ArrayList<>();
        ArrayList<Double> serviceStartTimeCurrentNode = new ArrayList<>();
        ArrayList<Double> serviceDurationCurrentNode = new ArrayList<>();
        ArrayList<Double> totalWaitTimeAfterLeavingCurrentNode = new ArrayList<>();
        ArrayList<Double> enterEDTime = new ArrayList<>();
        ArrayList<Double> leaveEDTime = new ArrayList<>();
        ArrayList<String> patientStatus = new ArrayList<>();
        ArrayList<Integer> patientConsultationDoctor = new ArrayList<>(); // 初始化为-1，表示尚未分配医生

        ArrayList<Double> firstTreatmentTime = new ArrayList<>();


        // 医生数量配置

        int numServerConRoom = alt[0];
        int numServerConCubicle = alt[1];
        int numServerResuscitation = alt[2];
        int numServerLabtest = alt[3];
        int numServerRegistration = alt[4];
        int numServerTriage = alt[5];

        // 为各个节点创建队列
        ArrayList<ArrayList<Integer>> conRoomQueue = new ArrayList<>();
        for (int i = 0; i < numServerConRoom; i++) {
            conRoomQueue.add(new ArrayList<>());
        }
        ArrayList<ArrayList<Integer>> conCubicleQueue = new ArrayList<>();
        for (int i = 0; i < numServerConCubicle; i++) {
            conCubicleQueue.add(new ArrayList<>());
        }
        ArrayList<Integer> resuscitationQueue = new ArrayList<>();
        ArrayList<Integer> labtestQueue = new ArrayList<>();
        ArrayList<Integer> registrationQueue = new ArrayList<>();
        ArrayList<Integer> triageQueue = new ArrayList<>();



        Map<String, ArrayList<Integer>> nodeQueueMap = new HashMap<>();
        for (int i = 0; i < conRoomQueue.size(); i++) {
            nodeQueueMap.put("conRoom" + (i), conRoomQueue.get(i));  ////////////////////
        }
        for (int i = 0; i < conCubicleQueue.size(); i++) {
            nodeQueueMap.put("conCubicle" + (i), conCubicleQueue.get(i));  ////////////////////
        }

        nodeQueueMap.put("resuscitation", resuscitationQueue);
        nodeQueueMap.put("labtest", labtestQueue);
        nodeQueueMap.put("registration", registrationQueue);
        nodeQueueMap.put("triage", triageQueue);


        Map<String, Integer> serverCount = new HashMap<>();
        serverCount.put("conRoom", numServerConRoom);
        serverCount.put("conCubicle", numServerConCubicle);
        serverCount.put("resuscitation", numServerResuscitation);
        serverCount.put("labtest", numServerLabtest);
        serverCount.put("registration", numServerRegistration);
        serverCount.put("triage", numServerTriage);

        // 仿真50名患者
        int numberOfPatients = 1500;
        double lambda = (double) 1 ;

        double currentTime = 0;
        ArrayList<Integer> patientInED = new ArrayList<>();

        for (int i = 0; i < numberOfPatients; i++) {
            patientId.add(i);
            double rand = R.nextDouble();
            if (rand < 0.01) {
                patientType.add("C1");
            }else if (rand < 0.03) {
                patientType.add("C2");
            }else if (rand < 0.40) {
                patientType.add("C3");
            }else if (rand < 0.95) {
                patientType.add("C4");
            }else {
                patientType.add("C5");
            }

            currentNode.add("outsideED");
            nextNode.add("enterED");

            double interArrivalTime = -Math.log(1.0 - R.nextDouble()) / lambda;
            double arrivalTime = currentTime + interArrivalTime;
            nextNodeTime.add(arrivalTime);
            enterEDTime.add(arrivalTime);

            // 更新当前时间
            currentTime += interArrivalTime;

            arrivalTimeCurrentNode.add(null);
            serviceStartTimeCurrentNode.add(null);
            serviceDurationCurrentNode.add(null);
            totalWaitTimeAfterLeavingCurrentNode.add(0.0);
            leaveEDTime.add(null);
            patientStatus.add("notEnteredED");
            patientConsultationDoctor.add(-1);
            firstTreatmentTime.add(null);
        }
        int count = 0;
        long ____t1 = System.nanoTime();

        while (count==0 || !patientInED.isEmpty()) {
            long ____t2 = System.nanoTime();
            int nextPatientIndex = -1;
            double earliestNextNodeTime = Double.MAX_VALUE;
            if(count == 0 ){
                patientInED.add(0);
                patientInED.add(1);
                earliestNextNodeTime = nextNodeTime.get(0);
                nextPatientIndex = 0;
                count = 1;
            }else{
                for(int i = 0; i < patientInED.size();i++){
                    if(nextNodeTime.get(patientInED.get(i)) < earliestNextNodeTime) {
                        earliestNextNodeTime = nextNodeTime.get(patientInED.get(i));
                        nextPatientIndex = patientInED.get(i);
                    }
                }
                if(nextPatientIndex == count && count < numberOfPatients-1){
                    patientInED.add(count+1);
                    count++;
                }
            }


//            if (nextPatientIndex == -1) {
//                System.out.println("12312312k3jl1k2j3lk12j3lkj1l2kj3l1k2jl3kj12lk3jl1");
//                break;
//            }


            currentTime = earliestNextNodeTime;

            //更新患者状态和节点
            if (patientStatus.get(nextPatientIndex).equals("notEnteredED")&& !Objects.equals(patientType.get(nextPatientIndex), "C1")) {
                currentNode.set(nextPatientIndex,"registration");
                nextNode.set(nextPatientIndex,"triage");
                arrivalTimeCurrentNode.set(nextPatientIndex,currentTime);
                patientStatus.set(nextPatientIndex, "inED");
                serviceDurationCurrentNode.set(nextPatientIndex,getProcessingTime(patientType.get(nextPatientIndex), "registration")); // 使用新的处理时间方法

                if (registrationQueue.size() < numServerRegistration) {
                    registrationQueue.add(nextPatientIndex);
                    serviceStartTimeCurrentNode.set(nextPatientIndex, currentTime);
                    nextNodeTime.set(nextPatientIndex, currentTime + serviceDurationCurrentNode.get(nextPatientIndex));
                } else {
                    double maxPreviousDepartureTime = getMaxPreviousDepartureTime(registrationQueue, registrationQueue.size(), numServerRegistration, nextNodeTime);
                    serviceStartTimeCurrentNode.set(nextPatientIndex, maxPreviousDepartureTime);
                    nextNodeTime.set(nextPatientIndex, maxPreviousDepartureTime + serviceDurationCurrentNode.get(nextPatientIndex));
                    registrationQueue.add(nextPatientIndex);
                }
            }else if(patientStatus.get(nextPatientIndex).equals("notEnteredED") && Objects.equals(patientType.get(nextPatientIndex), "C1")){
                currentNode.set(nextPatientIndex,"resuscitation");
                nextNode.set(nextPatientIndex,"labtest");
                arrivalTimeCurrentNode.set(nextPatientIndex,currentTime);
                patientStatus.set(nextPatientIndex, "inED");
                serviceDurationCurrentNode.set(nextPatientIndex,getProcessingTime(patientType.get(nextPatientIndex), "resuscitation")); // 使用新的处理时间方法
                addToQueue(nextPatientIndex, "resuscitation", nodeQueueMap, patientType, serverCount, nextNodeTime,
                        currentTime, serviceStartTimeCurrentNode, serviceDurationCurrentNode);

            }else if (nextNode.get(nextPatientIndex).equals("triage")){
                totalWaitTimeAfterLeavingCurrentNode.set(nextPatientIndex, totalWaitTimeAfterLeavingCurrentNode.get(nextPatientIndex)
                        + serviceStartTimeCurrentNode.get(nextPatientIndex) - arrivalTimeCurrentNode.get(nextPatientIndex));
                registrationQueue.remove(Integer.valueOf(nextPatientIndex));
                currentNode.set(nextPatientIndex, "triage");
                nextNode.set(nextPatientIndex, determineNextNode(patientType.get(nextPatientIndex),"triage"));
                arrivalTimeCurrentNode.set(nextPatientIndex, currentTime);
                serviceDurationCurrentNode.set(nextPatientIndex, getProcessingTime(patientType.get(nextPatientIndex), "triage"));

                if (triageQueue.size() < numServerTriage) {
                    triageQueue.add(nextPatientIndex);
                    serviceStartTimeCurrentNode.set(nextPatientIndex, currentTime);
                    nextNodeTime.set(nextPatientIndex, currentTime + serviceDurationCurrentNode.get(nextPatientIndex));
                } else {
                    double maxPreviousDepartureTime = getMaxPreviousDepartureTime(triageQueue, triageQueue.size(), numServerTriage, nextNodeTime);
                    serviceStartTimeCurrentNode.set(nextPatientIndex, maxPreviousDepartureTime);
                    nextNodeTime.set(nextPatientIndex, maxPreviousDepartureTime + serviceDurationCurrentNode.get(nextPatientIndex));
                    triageQueue.add(nextPatientIndex);
                }

            }else if (nextNode.get(nextPatientIndex).equals("discharge") && !Objects.equals(patientType.get(nextPatientIndex), "C1")){
                if (firstTreatmentTime.get(nextPatientIndex)==null){
                    firstTreatmentTime.set(nextPatientIndex,serviceStartTimeCurrentNode.get(nextPatientIndex)-enterEDTime.get(nextPatientIndex));
//                    System.out.println("2-ID: " + patientId.get(nextPatientIndex)+" T: " + patientType.get(nextPatientIndex)+" 设置了FT"+firstTreatmentTime.get(nextPatientIndex));
                }
                totalWaitTimeAfterLeavingCurrentNode.set(nextPatientIndex, totalWaitTimeAfterLeavingCurrentNode.get(nextPatientIndex)
                        + serviceStartTimeCurrentNode.get(nextPatientIndex) - arrivalTimeCurrentNode.get(nextPatientIndex));
                deleteFromQueueOfConsultation(nextPatientIndex, nodeQueueMap, patientConsultationDoctor, currentNode.get(nextPatientIndex));
                currentNode.set(nextPatientIndex,"discharge");
                nextNode.set(nextPatientIndex,"discharge");
                arrivalTimeCurrentNode.set(nextPatientIndex,currentTime);
                nextNodeTime.set(nextPatientIndex, Double.MAX_VALUE);
                patientStatus.set(nextPatientIndex, "leftED");

                patientInED.remove(Integer.valueOf(nextPatientIndex));//zhong

                serviceStartTimeCurrentNode.set(nextPatientIndex,currentTime);
                serviceDurationCurrentNode.set(nextPatientIndex, Double.MAX_VALUE);
                leaveEDTime.set(nextPatientIndex, currentTime);

            }else if (nextNode.get(nextPatientIndex).equals("discharge") && Objects.equals(patientType.get(nextPatientIndex), "C1")){

                totalWaitTimeAfterLeavingCurrentNode.set(nextPatientIndex, totalWaitTimeAfterLeavingCurrentNode.get(nextPatientIndex)
                        + serviceStartTimeCurrentNode.get(nextPatientIndex) - arrivalTimeCurrentNode.get(nextPatientIndex));
                deleteFromQueue(nextPatientIndex, currentNode.get(nextPatientIndex), nodeQueueMap);
                currentNode.set(nextPatientIndex,"discharge");
                nextNode.set(nextPatientIndex,"discharge");
                arrivalTimeCurrentNode.set(nextPatientIndex,currentTime);
                nextNodeTime.set(nextPatientIndex, Double.MAX_VALUE);
                patientStatus.set(nextPatientIndex, "leftED");

                patientInED.remove(Integer.valueOf(nextPatientIndex));//zhong

                serviceStartTimeCurrentNode.set(nextPatientIndex,currentTime);
                serviceDurationCurrentNode.set(nextPatientIndex, Double.MAX_VALUE);
                leaveEDTime.set(nextPatientIndex, currentTime);

            }else if (nextNode.get(nextPatientIndex).equals("conRoom") || nextNode.get(nextPatientIndex).equals("conCubicle")){

                if (Objects.equals(currentNode.get(nextPatientIndex), "resuscitation")){
                    if (firstTreatmentTime.get(nextPatientIndex)!=null){
                        System.out.println("计算去过急诊室的c1c2患者的ft时出现问题！！！");
                    }
                    firstTreatmentTime.set(nextPatientIndex,serviceStartTimeCurrentNode.get(nextPatientIndex)-enterEDTime.get(nextPatientIndex));
                }


                totalWaitTimeAfterLeavingCurrentNode.set(nextPatientIndex, totalWaitTimeAfterLeavingCurrentNode.get(nextPatientIndex)
                        + serviceStartTimeCurrentNode.get(nextPatientIndex) - arrivalTimeCurrentNode.get(nextPatientIndex));
                deleteFromQueue(nextPatientIndex, currentNode.get(nextPatientIndex), nodeQueueMap);
                currentNode.set(nextPatientIndex, nextNode.get(nextPatientIndex));
                nextNode.set(nextPatientIndex, determineNextNode(patientType.get(nextPatientIndex),currentNode.get(nextPatientIndex)));
                arrivalTimeCurrentNode.set(nextPatientIndex,currentTime);
                serviceDurationCurrentNode.set(nextPatientIndex,getProcessingTime(patientType.get(nextPatientIndex), currentNode.get(nextPatientIndex)));
                addToQueueOfConsultation(nextPatientIndex, nodeQueueMap, patientType, serverCount, nextNodeTime, currentTime, serviceStartTimeCurrentNode,
                        serviceDurationCurrentNode, patientConsultationDoctor, currentNode.get(nextPatientIndex));

            }else if (nextNode.get(nextPatientIndex).equals("labtest") && !Objects.equals(patientType.get(nextPatientIndex), "C1")){

                if (firstTreatmentTime.get(nextPatientIndex)==null){
                    firstTreatmentTime.set(nextPatientIndex,serviceStartTimeCurrentNode.get(nextPatientIndex)-enterEDTime.get(nextPatientIndex));
//                    System.out.println("5-ID: " + patientId.get(nextPatientIndex)+" T: " + patientType.get(nextPatientIndex)+" 设置了FT"+firstTreatmentTime.get(nextPatientIndex));
                }
                totalWaitTimeAfterLeavingCurrentNode.set(nextPatientIndex, totalWaitTimeAfterLeavingCurrentNode.get(nextPatientIndex)
                        + serviceStartTimeCurrentNode.get(nextPatientIndex) - arrivalTimeCurrentNode.get(nextPatientIndex));
                deleteFromQueueOfConsultation(nextPatientIndex, nodeQueueMap, patientConsultationDoctor, currentNode.get(nextPatientIndex));
                if (Objects.equals(currentNode.get(nextPatientIndex), "conRoom")){
                    nextNode.set(nextPatientIndex, "conRoom");
                }else if (Objects.equals(currentNode.get(nextPatientIndex), "conCubicle")){
                    nextNode.set(nextPatientIndex, "conCubicle");
                }
                currentNode.set(nextPatientIndex, "labtest");
                arrivalTimeCurrentNode.set(nextPatientIndex, currentTime);
                serviceDurationCurrentNode.set(nextPatientIndex, getProcessingTime(patientType.get(nextPatientIndex), currentNode.get(nextPatientIndex)));
                addToQueue(nextPatientIndex, currentNode.get(nextPatientIndex), nodeQueueMap, patientType, serverCount, nextNodeTime,
                        currentTime, serviceStartTimeCurrentNode, serviceDurationCurrentNode);

            }else{
                if (currentNode.get(nextPatientIndex)== "resuscitation"){
                    if (firstTreatmentTime.get(nextPatientIndex)!=null){
                        System.out.println("计算去过急诊室的c1c2患者的ft时出现问题！！！");
                    }
                    firstTreatmentTime.set(nextPatientIndex,serviceStartTimeCurrentNode.get(nextPatientIndex)-enterEDTime.get(nextPatientIndex));
                }
                totalWaitTimeAfterLeavingCurrentNode.set(nextPatientIndex, totalWaitTimeAfterLeavingCurrentNode.get(nextPatientIndex)
                        + serviceStartTimeCurrentNode.get(nextPatientIndex) - arrivalTimeCurrentNode.get(nextPatientIndex));
                deleteFromQueue(nextPatientIndex, currentNode.get(nextPatientIndex), nodeQueueMap);
                currentNode.set(nextPatientIndex, nextNode.get(nextPatientIndex));
                nextNode.set(nextPatientIndex, determineNextNode(patientType.get(nextPatientIndex),currentNode.get(nextPatientIndex)));
                arrivalTimeCurrentNode.set(nextPatientIndex, currentTime);
                serviceDurationCurrentNode.set(nextPatientIndex, getProcessingTime(patientType.get(nextPatientIndex), currentNode.get(nextPatientIndex)));
                addToQueue(nextPatientIndex, currentNode.get(nextPatientIndex), nodeQueueMap, patientType, serverCount, nextNodeTime,
                        currentTime, serviceStartTimeCurrentNode, serviceDurationCurrentNode);
            }
            long ____t3 = System.nanoTime();
            sortingTime = sortingTime - ____t3+____t2;
        }
        long ____t4 = System.nanoTime();
        sortingTime = sortingTime + ____t4 - ____t1;

        long endTime = System.currentTimeMillis();

        int numberOfC1 = 0;
        int numberOfC2 = 0;
        int numberOfC3 = 0;
        int numberOfC4 = 0;
        int numberOfC5 = 0;

        double totalWaitTimeOfC4 = 0.0;
        double totalWaitTimeOfC5 = 0.0;

        int numberOfGoodC1FT = 0;
        int numberOfGoodC2FT = 0;
        int numberOfGoodC3FT = 0;

        double goodC1FT = 1;
        double goodC2FT = 15;
        double goodC3FT = 30;

        double perOfSatisfyC1 = 0.99;
        double perOfSatisfyC2 = 0.95;
        double perOfSatisfyC3 = 0.9;

        int indexOfStartCounting = (int)(0.1*numberOfPatients);
        int indexOfEndCounting = (int)(0.9*numberOfPatients);
        for (int i = indexOfStartCounting; i < indexOfEndCounting; i++) {
            if (patientType.get(i).equals("C1")){
                numberOfC1 = numberOfC1 + 1;
                if (firstTreatmentTime.get(i)<= goodC1FT){
                    numberOfGoodC1FT = numberOfGoodC1FT + 1;
                }
            }else if(patientType.get(i).equals("C2")){
                numberOfC2 = numberOfC2 + 1;
                if (firstTreatmentTime.get(i)<= goodC2FT){
                    numberOfGoodC2FT = numberOfGoodC2FT + 1;
                }
            }else if(patientType.get(i).equals("C3")){
                numberOfC3 = numberOfC3 + 1;
                if (firstTreatmentTime.get(i)<= goodC3FT){
                    numberOfGoodC3FT = numberOfGoodC3FT + 1;
                }
            }else if(patientType.get(i).equals("C4")){
                numberOfC4 = numberOfC4 + 1;
                totalWaitTimeOfC4 = totalWaitTimeOfC4 + totalWaitTimeAfterLeavingCurrentNode.get(i);
            }else{
                numberOfC5 = numberOfC5 + 1;
                totalWaitTimeOfC5 = totalWaitTimeOfC5 + totalWaitTimeAfterLeavingCurrentNode.get(i);
            }
        }


        //ED4
        if ((double)numberOfGoodC1FT/numberOfC1>=perOfSatisfyC1 || numberOfC1==0){
            satisfyC1OrNot= 1.0;
        }else{
            satisfyC1OrNot= 0.0;
        }
        if ((double)numberOfGoodC2FT/numberOfC2>=perOfSatisfyC2 || numberOfC2==0){
            satisfyC2OrNot= 1.0;
        }else{
            satisfyC2OrNot= 0.0;
        }
        if ((double)numberOfGoodC3FT/numberOfC3>=perOfSatisfyC3 || numberOfC3==0){
            satisfyC3OrNot= 1.0;
        }else{
            satisfyC3OrNot= 0.0;
        }
        averageWTForC4C5 = (double)(totalWaitTimeOfC4+totalWaitTimeOfC5)/(numberOfC4+numberOfC5);

    }
    // 确定下一个节点的方法
    public static String determineNextNode(String patientType, String currentNode) {
        switch (currentNode) {
            case "enterED":
                if (patientType.equals("C1")) {
                    return "resuscitation";
                }else {
                    return "registration";
                }
            case "registration":
                return "triage";
            case "triage":
                if (patientType.equals("C5") || patientType.equals("C4")) {
                    return "conRoom";
                }else if (patientType.equals("C3") ){
                    return "conCubicle";
                }else if (patientType.equals("C2") ){
                    if (R.nextDouble() < 0.05) {
                        return "resuscitation";
                    }else{
                        return "conCubicle";
                    }
                }else{
                    System.out.println("分诊处出现了不该出现的患者类型！！！");
                }
            case "resuscitation":
                if (patientType.equals("C2")) {
                    return "conCubicle";
                }else if (patientType.equals("C1")){
                    return "labtest";
                }else{
                    System.out.println("急救处出现了不该出现的患者类型！！！");
                }
            case "conRoom":
                if(patientType.equals("C4")){
                    if (R.nextDouble() < 0.15) {
                        return "labtest";
                    } else {
                        return "discharge";
                    }
                }else if(patientType.equals("C5")){
                    if (R.nextDouble() < 0.03) {
                        return "labtest";
                    } else {
                        return "discharge";
                    }
                }else{
                    System.out.println("conRoom处出现了不该出现的患者类型！！！");
                }
            case "conCubicle":
                if(patientType.equals("C2")){
                    if (R.nextDouble() < 0.1) {
                        return "labtest";
                    } else {
                        return "discharge";
                    }
                }else if(patientType.equals("C3")){
                    if (R.nextDouble() < 0.2) {
                        return "labtest";
                    } else {
                        return "discharge";
                    }
                }else{
//                    System.out.println("conCubicle处出现了不该出现的患者类型！！！");
                }
            case "labtest":
                if (patientType.equals("C5") || patientType.equals("C4")) {
                    return "conRoom";
                }else if (patientType.equals("C3") || (patientType.equals("C2"))){
                    return "conCubicle";
                }else if (patientType.equals("C1")){
                    return "discharge";
                }else{
//                    System.out.println("labtest处出现了不该出现的患者类型！！！");
                }
            default:
//                System.out.println("判断下一地点时出现不正确的地点!!!");
                return "Unknown";
        }
    }

    public static double getProcessingTime(String patientType, String currentNode) {
        switch (currentNode) {
            case "enterED":
                return 0;
            case "registration":
                return triangularDistribution(1.5, 2, 2.5);
            case "triage":
                if (patientType.equals("C2")) {
                    return 3.5;
                }else if(patientType.equals("C3") || patientType.equals("C4")){
                    return triangularDistribution(2, 2.5, 3);
                }else if(patientType.equals("C5")){
                    return triangularDistribution(2, 2.25, 2.5);
                } else {
//                    System.out.println("time: triage处出现了不该出现的患者类型！！！");
                }
            case "resuscitation":
                if (patientType.equals("C1")) {
                    return triangularDistribution(30, 35, 45);
                }else if(patientType.equals("C2")){
                    return triangularDistribution(20, 30, 35);
                }else {
//                    System.out.println("time: resuscitation处出现了不该出现的患者类型！！！");
                }
            case "conRoom":
                if (patientType.equals("C5")) {
                    return triangularDistribution(8, 8.5, 9);
                }else if(patientType.equals("C4")){
                    return triangularDistribution(8.5, 9, 9.5);
                }else {
                    System.out.println("time: conRoom处出现了不该出现的患者类型！！！");
                }
            case "conCubicle":
                if (patientType.equals("C3")) {
                    return triangularDistribution(9, 9.5, 10);
                }else if(patientType.equals("C2")){
                    return triangularDistribution(9.5, 10, 10.5);
                }else {
//                    System.out.println("time: conCubicle处出现了不该出现的患者类型！！！");
                }
            case "labtest":
                if (patientType.equals("C1")) {
                    return triangularDistribution(0, 4, 8);
                }else if(patientType.equals("C2")){
                    return triangularDistribution(0, 2, 25);
                }else if(patientType.equals("C3")){
                    return triangularDistribution(0, 2, 45);
                }else if(patientType.equals("C4")){
                    return triangularDistribution(0, 2.5, 60);
                }else if(patientType.equals("C5")){
                    return triangularDistribution(2, 6, 60);
                }else {
//                    System.out.println("time: labtest处出现了不该出现的患者类型！！！");
                }
            default:
//                System.out.println("判断time时出现不正确的地点!!!");
                return 0;
        }
    }


    private static double triangularDistribution(double min, double mode, double max) {
        double f = (mode - min) / (max - min);
        double rand = R.nextDouble();

        if (rand < f) {
            return min + Math.sqrt(rand * (max - min) * (mode - min));
        } else {
            return max - Math.sqrt((1 - rand) * (max - min) * (max - mode));
        }
    }

    public static double uniformDistribution(double min, double max) {
        return min + (max - min) * R.nextDouble();
    }

    private static double exponentialDistribution(double mean) {
        return -mean * Math.log(1.0 - R.nextDouble());
    }

    public static void deleteFromQueue(int patientId, String currentNode, Map<String, ArrayList<Integer>> nodeQueueMap){
        ArrayList<Integer> lastQueue = nodeQueueMap.get(currentNode);
        lastQueue.remove(Integer.valueOf(patientId));
    }

    public static void deleteFromQueueOfConsultation(int patientId, Map<String, ArrayList<Integer>> nodeQueueMap,
                                                     ArrayList<Integer> patientConsultationDoctor, String currentNode){

        int doctorIndex;
        doctorIndex = patientConsultationDoctor.get(patientId);
        ArrayList<Integer> lastQueue = nodeQueueMap.get(currentNode+doctorIndex);
        lastQueue.remove(Integer.valueOf(patientId));
    }

    public void addToQueue(int patientId, String currentNode, Map<String, ArrayList<Integer>> nodeQueueMap,
                                  ArrayList<String> patientType, Map<String, Integer> serverCount, ArrayList<Double> nextNodeTime,
                                  double currentTime, ArrayList<Double> serviceStartTimeCurrentNode,
                                  ArrayList<Double> serviceDurationCurrentNode) {

        ArrayList<Integer> queue = nodeQueueMap.get(currentNode);
        int doctorsAvailable = serverCount.get(currentNode);

        String type = patientType.get(patientId);
        int priority = getTypePriority(type);

        if (queue.size() < doctorsAvailable) {
            queue.add(patientId);
            serviceStartTimeCurrentNode.set(patientId, currentTime);
            nextNodeTime.set(patientId, currentTime + serviceDurationCurrentNode.get(patientId));
        } else {
            int insertIndex = queue.size();
            for (int i = doctorsAvailable; i < queue.size(); i++) {
                int existingPatientId = queue.get(i);
                int existingPatientPriority = getTypePriority(patientType.get(existingPatientId));

                if (priority < existingPatientPriority) {
                    insertIndex = i;
                    break;
                }
            }
            queue.add(insertIndex, patientId);

            for (int j = insertIndex; j < queue.size(); j++) {
                int currentPatientId = queue.get(j);
//                long ____t = System.nanoTime();
                double maxPreviousDepartureTime = getMaxPreviousDepartureTime(queue, j, doctorsAvailable, nextNodeTime);
//                sortingTime = sortingTime + System.nanoTime() - ____t;
                serviceStartTimeCurrentNode.set(currentPatientId, maxPreviousDepartureTime);
                nextNodeTime.set(currentPatientId, maxPreviousDepartureTime + serviceDurationCurrentNode.get(currentPatientId));
            }
//            for (int i = 0; i < queue.size(); i++){
//                System.out.println("                "+i+" "+queue.get(i)+" "+patientType.get(queue.get(i))+" "+nextNodeTime.get(queue.get(i)));
//            }

        }
    }

    public static int getDoctorWithMinQueue(String currentNode, Map<String, ArrayList<Integer>> nodeQueueMap, int numDoctors) {
        int minQueueSize = Integer.MAX_VALUE;
        int selectedDoctor = -1;
        for (int i = 0; i < numDoctors; i++) {
            int queueSize = nodeQueueMap.get(currentNode + i).size(); //????????
            if (queueSize < minQueueSize) {
                minQueueSize = queueSize;
                selectedDoctor = i;
            }
        }
        return selectedDoctor;
    }

    public static void addToQueueOfConsultation(int patientId, Map<String, ArrayList<Integer>> nodeQueueMap,
                                                ArrayList<String> patientType, Map<String, Integer> serverCount, ArrayList<Double> nextNodeTime,
                                                double currentTime, ArrayList<Double> serviceStartTimeCurrentNode,
                                                ArrayList<Double> serviceDurationCurrentNode,
                                                ArrayList<Integer> patientConsultationDoctor, String currentNode) {
        int doctorIndex = -1;
        if (patientConsultationDoctor.get(patientId) != -1) {
            doctorIndex = patientConsultationDoctor.get(patientId);
        } else {
            if (currentNode.equals("conRoom")){
                doctorIndex = getDoctorWithMinQueue(currentNode, nodeQueueMap, serverCount.get("conRoom"));
                patientConsultationDoctor.set(patientId, doctorIndex);   //doctorindex 0,1,2,3,...,x
            }else if(currentNode.equals("conCubicle")){
                doctorIndex = getDoctorWithMinQueue(currentNode, nodeQueueMap, serverCount.get("conCubicle"));
                patientConsultationDoctor.set(patientId, doctorIndex);   //doctorindex 0,1,2,3,...,x
            }
        }
        ArrayList<Integer> queue = nodeQueueMap.get(currentNode+doctorIndex);

        String type = patientType.get(patientId);
        int priority = getTypePriority(type);

        if (queue.isEmpty()) {
            queue.add(patientId);
            serviceStartTimeCurrentNode.set(patientId, currentTime);
            nextNodeTime.set(patientId, currentTime + serviceDurationCurrentNode.get(patientId));
        } else {

            int insertIndex = queue.size();
            for (int i = 1; i < queue.size(); i++) {
                int existingPatientId = queue.get(i);
                int existingPatientPriority = getTypePriority(patientType.get(existingPatientId));

                if (priority < existingPatientPriority) {
                    insertIndex = i;
                    break;
                }
            }
            queue.add(insertIndex, patientId);

            for (int j = insertIndex; j < queue.size(); j++) {
                int currentPatientId = queue.get(j);
                double maxPreviousDepartureTime = nextNodeTime.get(queue.get(j - 1));
                serviceStartTimeCurrentNode.set(currentPatientId , maxPreviousDepartureTime);
                nextNodeTime.set(currentPatientId, maxPreviousDepartureTime + serviceDurationCurrentNode.get(currentPatientId));
            }
        }
    }

    public static double getMaxPreviousDepartureTime(ArrayList<Integer> queue, int currentIndex, int doctorsAvailable, ArrayList<Double> nextNodeTime) {

        if (doctorsAvailable <= 0) {
            throw new IllegalArgumentException("doctorsAvailable must be >= 1");
        }

        PriorityQueue<Double> minHeap = new PriorityQueue<>(doctorsAvailable);

        for (int i = 0; i < currentIndex; i++) {
            double time = nextNodeTime.get(queue.get(i));

            if (minHeap.size() < doctorsAvailable) {
                minHeap.offer(time);
            } else if (time > minHeap.peek()) {
                minHeap.poll();
                minHeap.offer(time);
            }
        }

        if (minHeap.size() < doctorsAvailable) {
            System.out.println("Warning: Less than doctorsAvailable elements in history!");
        }

        return minHeap.peek();  // 第 k 大的元素（即最小堆的堆顶）


//        if (doctorsAvailable <= 0) {
//            throw new IllegalArgumentException("doctorsAvailable must be >= 1");
//        }
//
//        PriorityQueue<Double> minHeap = new PriorityQueue<>(doctorsAvailable);
//
//        for (int i = 0; i < currentIndex; i++) {
//            double time = nextNodeTime.get(queue.get(i));
//
//            if (minHeap.size() < doctorsAvailable) {
//                minHeap.offer(time);
//            } else if (time > minHeap.peek()) {
//                minHeap.poll();
//                minHeap.offer(time);
//            }
//        }
//
//        if (minHeap.size() < doctorsAvailable) {
//            System.out.println("Warning: Less than doctorsAvailable elements in history!");
//        }
//
//        return minHeap.peek();  // 第 k 大的元素（即最小堆的堆顶）


    }


    // 获取患者类型的优先级
    public static int getTypePriority(String type) {
        switch (type) {
            case "C1":
                return 1;
            case "C2":
                return 2;
            case "C3":
                return 3;
            case "C4":
                return 4;
            case "C5":
                return 5;
            default:
//                System.out.println("判断患者优先级时出现问题！！！");
                return Integer.MAX_VALUE;

        }
    }
}




