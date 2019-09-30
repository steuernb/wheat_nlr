
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import support.BioSequence;
import support.FastaReader;

/**
 * 
 * These are the scripts for creating the figure (currently Figure 3) that shows how deep you have to sequence to see a certain amount of NLRs.
 * 
 * 
 * @author steuernb
 *
 */
public class NLR_Expression_Level {

	
	public static void main(String[] args){
		try {
			
			
			File nlrTranscriptsFile = new File("steuernb/nlr_paper/v1.1_analysis/intermediate/NLR-transcripts.txt");
			File hcCDS = new File("steuernb/nlr_paper/v1.1_analysis/inputData/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_cds.fasta");
			File lcCDS = new File("steuernb/nlr_paper/v1.1_analysis/inputData/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_LC_20170706_cds.fasta");
			File transcriptFile = new File("steuernb/nlr_paper/v1.1_analysis/rnaseq/NLR/ByGene/cDNA_RenSeq_count.tsv");
			File outputFile = new File("steuernb/nlr_paper/v1.1_analysis/intermediate/ExpressionPerGb.txt");
			
			HashMap<String, Integer> nlrs = readNLRTranscriptList(nlrTranscriptsFile, hcCDS, lcCDS);
			
			HashMap<String, double[]> counts = readFile(transcriptFile,  nlrs);
			
			BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
			out.write("modificator\tnumNLRs\tnumExpressed\tnumRNA\tnumRenSeq\tnumTotalRNA");
			out.newLine();
			printExpressedNLRs( counts) ;
			
			for( int i = 306; i>0; i--){
				runIteration(counts,i, out);
			}
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
	}
	
	
	public static void printExpressedNLRs( HashMap<String, double[]> counts ) {
		
		
		int detected = 0;
		int detectedRNA = 0;
		int detectedRNAOnly = 0;
		int detectedRenSeq = 0;
		int detectedRenSeqOnly = 0;
		
		for(Iterator<String> iterator = counts.keySet().iterator(); iterator.hasNext();) {
			
			String gene =iterator.next();
			double[] a = counts.get(gene);
			
			double level = 5*a[0];
			
			double e1 = (a[4] + a[5] + a[6]) *150;    //expression RNA  : sum of counts in all replicates times read length
			double e2 = (a[1] + a[2] + a[3]) *150;    // expression RenSeq
			
			if( e1 >=level|| e2 >= level) { detected++;}
			if( e1 >= level) {detectedRNA++;}	
			if( e2 >= level) {detectedRenSeq++;}
			
			if(e1>=level && e2 < level) {detectedRNAOnly++;}
			if(e2>=level && e1 < level) {detectedRenSeqOnly++;}
			
		}
		
		System.out.println("Number of detected genes (RNA/RenSeq): " + detected + " (" + detectedRNA + "/" + detectedRenSeq+")" );
		System.out.println("Number of genes detected in only one dataset (RNA/RenSeq): " +detectedRNAOnly+"/"+  detectedRenSeqOnly);
	}
	
	public static void runIteration(HashMap<String, double[]> h, int modificator, BufferedWriter out)throws IOException {
		
			
				
			//int modificator = 200; // /3.05559E+11 * 1000000000 * modificator
			
			// total rna-seq data that was used = 3.05559E+11 bp (~305 Gb)
			
			
			
			int countAll = 0;
			int countExpressed = 0;
			int countExpressedRNA = 0;
			int countRenSeq = 0;
			int countExpressedTotalRNA=0;
			HashSet<String> expressedGenes = new HashSet<String>();
			for(Iterator<String> iterator = h.keySet().iterator(); iterator.hasNext();){
				
				String id = iterator.next();
				double[] d = h.get(id);
				int length = (int) Math.round(d[0]);
				double minExpression = length * 5 / 150;
				double renseq = d[1] + d[2] + d[3];
				double rnaseq = d[4] + d[5] + d[6];
				
				if(rnaseq >=minExpression) {
					countExpressedTotalRNA++;
				}
				
				countAll ++;
				if( renseq >= minExpression || rnaseq >= minExpression){
					countExpressed++;
					
					expressedGenes.add(id.split("\\.")[0]);
					rnaseq = (rnaseq /  3.05559E+11) * 1000000000 * modificator;
					
					if( rnaseq >= minExpression){
						countExpressedRNA++;
					}
					if(renseq>=minExpression) {
						countRenSeq ++;
					}
					
					
				}
				
				
				
			}
			
			out.write(modificator + "\t" + countAll + "\t" + countExpressed  + "\t" + countExpressedRNA + "\t" + countRenSeq+"\t"+ countExpressedTotalRNA);
			out.newLine();
	}
	
	
	
	/**
	 * 
	 * read in the gene count file and the lengths of the longest transcripts. extract counts for NLRs for 3LE tissues in RenSeq and not enriched
	 * 
	 * @param inputFile
	 * @param nlrs
	 * @return
	 * @throws IOException
	 */
	public static HashMap<String, double[]> readFile(File inputFile, HashMap<String, Integer> nlrs)throws IOException{
		HashMap<String, double[]> h = new HashMap<String, double[]>();
		
		
		BufferedReader in = new BufferedReader(new FileReader(inputFile));

		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String[] split = inputline.split("\t");
			String gene = split[0];
			if( nlrs.containsKey(gene)){
				double[] d = new double[7];
				
				d[0] = nlrs.get(split[0]);
				d[1] = Double.parseDouble(split[1]);
				d[2] = Double.parseDouble(split[2]);
				d[3] = Double.parseDouble(split[3]);
				d[4] = Double.parseDouble(split[18]);
				d[5] = Double.parseDouble(split[19]);
				d[6] = Double.parseDouble(split[20]);
				h.put(split[0], d);
			}
		}

		in.close();
		
		
		return h ;
	}
	/**
	 * 
	 * Get one transcript per gene. Always the longest.
	 * 
	 * @param nlrTranscriptsFile
	 * @param inputCDS
	 * @return
	 * 		A Hashmap with genes as keys, length of longest transcript as value.
	 * @throws IOException
	 */
	private static HashMap<String, Integer> readNLRTranscriptList(File nlrTranscriptsFile, File hcCDS, File lcCDS)throws IOException{
		
		
		HashMap<String, Integer> h = new HashMap<String, Integer>();

		
		/*
		 * read in the names of all nlr transcripts 
		 */
		BufferedReader in = new BufferedReader(new FileReader(nlrTranscriptsFile));
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			h.put(inputline.split("\\.")[0].trim(), 0);  // reduce transcript name to gene name
		}
		in.close();
		
		
		/*
		 *  read in the sequences of all transcripts. Keep the length of the longest transcript per nlr gene
		 */
		File[] cds = {hcCDS, lcCDS};
		for( int i = 0; i< cds.length; i++) {
			FastaReader fastaReader = new FastaReader(cds[i]);
			for (BioSequence seq = fastaReader.readEntry(); seq != null; seq = fastaReader.readEntry()) {
				String gene = seq.getIdentifier().split("\\.")[0];
				
				if( h.containsKey(gene)) {
					if( h.get(gene) < seq.getLength()) {
						h.put(gene, seq.getLength());
					}
				}
				
			}
			fastaReader.close();
		}
		
		return h;
	}
	
}
