package brandeWulff.csNLRpaper;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;

import support.BioSequence;
import support.FastaReader;

public class PAMP_RnaSeq {

	public static void main(String[] args) {
		try {
			File hcCDS = new File("steuernb/nlr_paper/v1.1_analysis/inputData/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_cds.fasta");
			File lcCDS = new File("steuernb/nlr_paper/v1.1_analysis/inputData/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_LC_20170706_cds.fasta");
			File nlrTranscriptsFile = new File("steuernb/nlr_paper/v1.1_analysis/intermediate/NLR-transcripts.txt");
			
			
			
			File rawCounts = new File("steuernb/nlr_paper/v1.1_analysis/rnaseq/NLR/ByGene/PAMP_Triggered_Imune_Response_count.tsv");
			File rawTPM = new File("steuernb/nlr_paper/v1.1_analysis/rnaseq/NLR/ByGene/PAMP_Triggered_Imune_Response_tpm.tsv");
			
			File nlrTPM = new File("steuernb/nlr_paper/v1.1_analysis/rnaseq/NLR/ByGene/PAMP_Triggered_Imune_Response_tpm_nlrAnnotated.tsv");
			File degustFile = new File("steuernb/nlr_paper/v1.1_analysis/rnaseq/dge/degust.csv");
			
			
			File outputF = new File("steuernb/nlr_paper/v1.1_analysis/rnaseq/dge/compare_F.txt");
			File outputC = new File("steuernb/nlr_paper/v1.1_analysis/rnaseq/dge/compare_C.txt");
			File outputH = new File("steuernb/nlr_paper/v1.1_analysis/rnaseq/dge/compare_H.txt");
			File outputHCF = new File("steuernb/nlr_paper/v1.1_analysis/rnaseq/dge/compare_HCF.txt");
			
			
			//HashSet<String> keepGenes = filterCountData(rawCounts, hcCDS, lcCDS);
			HashSet<String> nlrs = getNLRGenes(nlrTranscriptsFile);
			
			
			writeNlrTpm(rawTPM, nlrs,  nlrTPM);
			
			
			//writeNLRsForTriangles(degustFile, rawTPM, nlrs, outputF, outputC, outputH, outputHCF);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	
	
	public static void writeNLRsForTriangles(File degustFile,File rawTPM, HashSet<String> nlrs
			, File outputF
			, File outputC
			, File outputH
			, File outputHCF
			
			)throws IOException{
		
		int numNLRsRead=0;
		
		/*
		 * read in data
		 */
		
		BufferedReader in = new BufferedReader(new FileReader(degustFile));
		in.readLine();
		
		HashMap<String, double[]> cpms = new HashMap<String, double[]>();	
		
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String[] split = inputline.split(",");
			double fdr = Double.parseDouble(split[8]);
			
			if(nlrs.contains(split[0])) {
				numNLRsRead++;
			}else {
				continue;
			}
			
			
			if( fdr > 0.05) { // filter for NLR genes and throw out genes with fdr of >0.05
				continue;
			}
			
			//32 is the first cpm column. HCF
			int index = 32;
			double h000 = (Double.parseDouble(split[index]) + Double.parseDouble(split[index+1]) + Double.parseDouble(split[index+2]))/3; index = index + 3;
			double h030 = (Double.parseDouble(split[index]) + Double.parseDouble(split[index+1]) + Double.parseDouble(split[index+2]))/3; index = index + 3;
			double h180 = (Double.parseDouble(split[index]) + Double.parseDouble(split[index+1]) + Double.parseDouble(split[index+2]))/3; index = index + 3;
			double c030 = (Double.parseDouble(split[index]) + Double.parseDouble(split[index+1]) + Double.parseDouble(split[index+2]))/3; index = index + 3;
			double c180 = (Double.parseDouble(split[index]) + Double.parseDouble(split[index+1]) + Double.parseDouble(split[index+2]))/3; index = index + 3;
			double f030 = (Double.parseDouble(split[index]) + Double.parseDouble(split[index+1]) + Double.parseDouble(split[index+2]))/3; index = index + 3;
			double f180 = (Double.parseDouble(split[index]) + Double.parseDouble(split[index+1]) + Double.parseDouble(split[index+2]))/3; index = index + 3;
			
			double[] a = {h000,h030,h180,c030,c180,f030,f180};
			
			cpms.put(split[0], a);
			
		}

		in.close();
		
		HashMap<String, double[] > tpms = new HashMap<String, double[]>();
		
		in = new BufferedReader(new FileReader(rawTPM));
		in.readLine();
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			
			String[] split = inputline.split("\t");
			if( !cpms.containsKey(split[0])) {
				continue;
			}
			
			int index = 1;
			double h000 = (Double.parseDouble(split[index]) +  Double.parseDouble(split[index+1]) + Double.parseDouble(split[index+2]))/3;index = index + 3;
			double h030 = (Double.parseDouble(split[index]) +  Double.parseDouble(split[index+1]) + Double.parseDouble(split[index+2]))/3;index = index + 3;
			double h180 = (Double.parseDouble(split[index]) +  Double.parseDouble(split[index+1]) + Double.parseDouble(split[index+2]))/3;index = index + 3;
			double c030 = (Double.parseDouble(split[index]) +  Double.parseDouble(split[index+1]) + Double.parseDouble(split[index+2]))/3;index = index + 3;
			double c180 = (Double.parseDouble(split[index]) +  Double.parseDouble(split[index+1]) + Double.parseDouble(split[index+2]))/3;index = index + 3;
			double f030 = (Double.parseDouble(split[index]) +  Double.parseDouble(split[index+1]) + Double.parseDouble(split[index+2]))/3;index = index + 3;
			double f180 = (Double.parseDouble(split[index]) +  Double.parseDouble(split[index+1]) + Double.parseDouble(split[index+2]))/3;index = index + 3;
			double[] a = {h000,h030,h180,c030,c180,f030,f180};
			tpms.put(split[0], a);
			
		}

		in.close();
		
		
		
		/*
		 * write data for triangles
		 */
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputF));
		out.write("gene\tH000\tF030\tF180\tmaxTPM");
		out.newLine();
		for(Iterator<String> iterator = cpms.keySet().iterator(); iterator.hasNext();) {
			String gene = iterator.next();
			out.write(gene + "\t" + cpms.get(gene)[0] + "\t" + cpms.get(gene)[5] + "\t" + cpms.get(gene)[6]);  //write the gene and the CPM values
			out.write("\t" + ( Math.max(Math.max(tpms.get(gene)[0], tpms.get(gene)[5]), tpms.get(gene)[6])  ) ); // write the max TPM value for the samples
			out.newLine();
		}
		out.close();
		
		
		out = new BufferedWriter(new FileWriter(outputC));
		out.write("gene\tH000\tC030\tC180\tmaxTPM");
		out.newLine();
		for(Iterator<String> iterator = cpms.keySet().iterator(); iterator.hasNext();) {
			String gene = iterator.next();
			out.write(gene + "\t" + cpms.get(gene)[0] + "\t" + cpms.get(gene)[3] + "\t" + cpms.get(gene)[4]);  //write the gene and the CPM values
			out.write("\t" + ( Math.max(Math.max(tpms.get(gene)[0], tpms.get(gene)[3]), tpms.get(gene)[4])  ) ); // write the max TPM value for the samples
			out.newLine();
		}
		out.close();
		
		
		 out = new BufferedWriter(new FileWriter(outputH));
		 out.write("gene\tH000\tH030\tH180\tmaxTPM");
			out.newLine();
			
		 for(Iterator<String> iterator = cpms.keySet().iterator(); iterator.hasNext();) {
			String gene = iterator.next();
			out.write(gene + "\t" + cpms.get(gene)[0] + "\t" + cpms.get(gene)[1] + "\t" + cpms.get(gene)[2]);  //write the gene and the CPM values
			out.write("\t" + ( Math.max(Math.max(tpms.get(gene)[0], tpms.get(gene)[1]), tpms.get(gene)[2])  ) ); // write the max TPM value for the samples
			out.newLine();
		}
		out.close();
		
		out = new BufferedWriter(new FileWriter(outputHCF));
		out.write("gene\tH030\tC030\tF030\tmaxTPM");
		out.newLine();
		
		for(Iterator<String> iterator = cpms.keySet().iterator(); iterator.hasNext();) {
			String gene = iterator.next();
			out.write(gene + "\t" + cpms.get(gene)[1] + "\t" + cpms.get(gene)[3] + "\t" + cpms.get(gene)[5]);  //write the gene and the CPM values
			out.write("\t" + ( Math.max(Math.max(tpms.get(gene)[1], tpms.get(gene)[3]), tpms.get(gene)[5])  ) ); // write the max TPM value for the samples
			out.newLine();
		}
		out.close();
		
		
		/*
		 * get differentially expressed genes.
		 * 
		 */
		
		int numF=0;
		int numC=0;
		int numBoth=0;
		int numEither=0;
		
		for(Iterator<String> iterator= cpms.keySet().iterator(); iterator.hasNext();) {
			String gene = iterator.next();
			double[] a = cpms.get(gene);
			//h000,h030,h180,c030,c180,f030,f180
			// 0    1     2    3    4   5     6
			if(a[3] >= 2*a[0]) { numC++;}
			if(a[5] >= 2*a[0]) { numF++;}
			if(a[3] >= 2*a[0]  && a[5] >= 2*a[0] ) { numBoth++;}
			if(a[3] >= 2*a[0]  || a[5] >= 2*a[0] ) { numEither++;}
		}
		System.out.println("number of genes interrogated: " +numNLRsRead);
		//System.out.println("number of genes interrogated: " + cpms.size());
		System.out.println("number of deg in F: " + numF);
		System.out.println("number of deg in C: " + numC);
		System.out.println("number of deg in both: " + numBoth);
		System.out.println("number of deg in either: " + numEither);
		
	}
	
	
	
	
	
	public static void writeNlrTpm(File rawTPM, HashSet<String> nlrs, File nlrTPM)throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(nlrTPM));
		BufferedReader in = new BufferedReader(new FileReader(rawTPM));
		out.write(in.readLine() + "\tNLR"  );
		out.newLine();
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String[] split = inputline.split("\t");
			
			Vector<Double> v= new Vector<Double>();
			for( int i = 1 ;i<split.length; i++) {
				v.add(Double.parseDouble(split[i]));
			}
			
			int num05 = 0;
			double min=v.get(0);
			double sum = 0;
			double max = 0;
			for(Iterator<Double> iterator = v.iterator(); iterator.hasNext();) {
				double d = iterator.next();
				if(d >0.5) {num05++;}
				sum = sum + d;
				if(d<min) {min=d;}
				if(d>max) {max=d;}
				
			}
			
			if (num05==0){
				continue;
			}
			if( max/min <2) {
				continue;
			}
			if( sum/v.size() < 0.5) {
				continue;
			}
			
			if( nlrs.contains(split[0]) ) {
				out.write(inputline + "\tNLR");
				out.newLine();
			}else {
				out.write(inputline + "\tno");
				out.newLine();
		
			}
		}

		in.close();
		out.close();
		
	}
	
	
	
	/**
	 * 
	 * read the list of nlr transcripts. return a HashSet with the gene names. These are gene names, not transcripts. If there is a "." in the name, everything after will be deleted!
	 * 
	 * @param nlrTranscriptsFile
	 * @return
	 * @throws IOException
	 */
	public static HashSet<String> getNLRGenes(File nlrTranscriptsFile)throws IOException{
		HashSet<String> nlrs = new HashSet<String>();
		BufferedReader in = new BufferedReader(new FileReader(nlrTranscriptsFile));
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			nlrs.add(inputline.split("\\.")[0].trim());  // reduce transcript name to gene name
		}
		in.close();
		return nlrs;
	}
	
	
	
	/**
	 * 
	 * check if a gene is sufficiently expressed by having longest transcript 5x covered (sum_count_of_replicates *150 >= longest_transcript *5)
	 * 
	 * @param inputFile
	 * 	raw count file
	 * @param hcCDS
	 * 	high confidence CDS sequences
	 * @param lcCDS
	 * 	low confidence CDS sequences
	 * @return
	 * 	HashSet with gene names that are to be kept
	 * @throws IOException
	 */
	public static HashSet<String> filterCountData(File inputFile,  File hcCDS, File lcCDS)throws IOException{
		
		HashMap<String,Integer> lengths = new HashMap<String, Integer>();
		File[] files = {hcCDS, lcCDS};
		for( int i = 0; i< files.length; i++) {
			File file =files[i];
			FastaReader fastaReader = new FastaReader(file);
			for (BioSequence seq = fastaReader.readEntry(); seq != null; seq = fastaReader.readEntry()) {
				String gene = seq.getIdentifier().split("\\.")[0];
				if(!lengths.containsKey(gene) || lengths.get(gene) < seq.getLength()){
					lengths.put(gene, seq.getLength());
				}
			}
			fastaReader.close();
		}
		
		HashSet<String> keepGenes = new HashSet<String>();
		
		
		BufferedReader in = new BufferedReader(new FileReader(inputFile));
		in.readLine();
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String[] split = inputline.split("\t");
			boolean keep = false;
			
			for( int i = 1; i< split.length; i = i+3) {
				double count = Double.parseDouble(split[i]) + Double.parseDouble(split[i+1]) + Double.parseDouble(split[i+2]) ;
				
				/*
				 * gene is considered expressed if replicates together have more than 5x coverage over the transcript
				 */
				if( count * 150 >= lengths.get(split[0])*5 ) {  
					keep = true;
				}
				
			}
			if(keep) {
				keepGenes.add(split[0]);
			}
			
		}

		in.close();
		
		return keepGenes;
		
		
		
	}
	
	
	
	
}