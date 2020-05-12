//    public static void writeBashCode() throws IOException{
//
//        File outputFile = new File(bashScriptPath);
//        if (!outputFile.exists()) {
//            outputFile.mkdir();
//        }
//
//        PrintWriter writer = new PrintWriter(new File(bashScriptPath + bashScriptFileName));
//
//        writer.println("#!/bin/bash");
//
//        writer.println("pwd");
//
//        writer.close();
//    }
//
//    public static void runBashCode() throws IOException {
//
//        writeBashCode();
//
//        String command = "sh " + bashScriptPath + bashScriptFileName;
//        Process p = Runtime.getRuntime().exec(command);
//
//        String s = null;
//
//        BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
//        System.out.println("Terminal output:");
//        while((s=stdInput.readLine()) != null) {
//            System.out.println(s);
//        }
//
//        System.out.println();
//
//        BufferedReader stdError = new BufferedReader(new InputStreamReader(p.getErrorStream()));
//        System.out.println("Terminal error:");
//        while((s=stdError.readLine()) != null) {
//            System.out.println(s);
//        }
//
//        System.exit(0);
//    }
