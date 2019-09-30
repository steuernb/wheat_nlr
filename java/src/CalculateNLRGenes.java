package brandeWulff.csNLRpaper;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;

import support.BioSequence;
import support.FastaReader;

public class CalculateNLRGenes {

	public static void main(String[] args) {
		try {
			
			/*
			 * ******** Input Data ************
			 */
			
			// IWGSC gene annotation v. 1.1
			File hcGFF = new File("steuernb/nlr_paper/v1.1_analysis/inputData/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.gff3");
			File lcGFF = new File("steuernb/nlr_paper/v1.1_analysis/inputData/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_LC_20170706.gff3");
			File hcCDS = new File("steuernb/nlr_paper/v1.1_analysis/inputData/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_cds.fasta");
			File lcCDS = new File("steuernb/nlr_paper/v1.1_analysis/inputData/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_LC_20170706_cds.fasta");
			
			
			// NLR Annotator output for CS Pseudomolecules
			File nlrTxt = new File("steuernb/nlr_paper/v1.1_analysis/NLR-Annotator/IWGSC.nlr.txt");
			File nlrsFasta = new File("steuernb/nlr_paper/v1.1_analysis/NLR-Annotator/cs.nlr.fasta");
			
			// NLR Parser (v1.0) on CS peptides v. 1.1
			File nlrPeptides= new File("steuernb/nlr_paper/v1.1_analysis/pep-mast/hc-and-lc-pep.nlrparser.txt");
			
			
			//Others
			File manualInspection = new File("steuernb/nlr_paper/v1.1_analysis/inputData/manual_inspection_nlrloci.txt");
			File adr1 = new File("steuernb/nlr_paper/v1.1_analysis/inputData/ADR1_addon.txt");
			/* 
			 ************************************** 
			 */
				
			
			
			/*
			 * ******** Output Data ************
			 */
				File tableS4 = new File("steuernb/nlr_paper/v1.1_analysis/outputData/tableS4.txt");
				File tableS5 = new File("steuernb/nlr_paper/v1.1_analysis/outputData/tableS5.txt");
				File nlrList = new File("steuernb/nlr_paper/v1.1_analysis/outputData/nlrList.txt");
			/* 
			 ************************************** 
			 */
				
			
			Hashtable<String, Hashtable<String, int[]>> nlrs = getNLRs(nlrTxt);
			Hashtable<String, Hashtable<String, int[]>> genes = getGFF(hcGFF, lcGFF);
			
			Hashtable<String, String[]> nlrInfo = getNLRInfo(nlrTxt);
			
			Hashtable<String, Hashtable<String, String>> peptides = readNLRParserPeptides(nlrPeptides);
			
			writeTableS4(nlrs, genes, nlrInfo, peptides, nlrsFasta, manualInspection,tableS4);
			writeTableS5(peptides, genes, nlrs,adr1, tableS5);
			
			
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	
	
	
		
		
	
	
	public static void writeTableS5(
			Hashtable<String, Hashtable<String, String>> peptides,
			Hashtable<String, Hashtable<String, int[]>> genes,
			Hashtable<String, Hashtable<String, int[]>> nlrs,
			File adr1,
			File outputFile
			)throws IOException{
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		out.write("gene id\ttranscript id\toverlapping NLR loci\ttranscript position\thas p-loop\thas NB-ARC motifs\thas LRRs\tmotifs");
		out.newLine();
		
		
		for(Enumeration<String> myenum1 = genes.keys(); myenum1.hasMoreElements();) {
			String chromosome = myenum1.nextElement();
			
			for(Enumeration<String> myenum2 = genes.get(chromosome).keys(); myenum2.hasMoreElements();) {
				
				String transcript = myenum2.nextElement();
				String gene = transcript.split("\\.")[0];
				
				int transcriptStart = genes.get(chromosome).get(transcript)[0];
				int transcriptEnd = genes.get(chromosome).get(transcript)[1];
				int transcriptStrand = genes.get(chromosome).get(transcript)[2];
				String position = chromosome + ":" + transcriptStart + "-" + transcriptEnd;
 				
				String overlappingNLRs = "";
				
				for(Enumeration<String> myenum3 = nlrs.get(chromosome).keys(); myenum3.hasMoreElements();) {
					String nlr = myenum3.nextElement();
					int nlrStart = nlrs.get(chromosome).get(nlr)[0];
					int nlrEnd = nlrs.get(chromosome).get(nlr)[1];
					int nlrStrand = nlrs.get(chromosome).get(nlr)[2];
					
					if( nlrStrand == transcriptStrand && Math.max(transcriptEnd, nlrEnd) - Math.min(transcriptStart, nlrStart)  < transcriptEnd-transcriptStart + nlrEnd-nlrStart  ) {  //right strand and overlapping
						overlappingNLRs = overlappingNLRs + "," + nlr;
					
					}
				}
				
				if( overlappingNLRs.length()==0 && (!peptides.containsKey(gene)  ||  !peptides.get(gene).containsKey(transcript) ) ) {
					continue;
				}
				
				
				if(overlappingNLRs.length() >0) {
					overlappingNLRs = overlappingNLRs.substring(1);
				}
				
				out.write(gene + "\t" + transcript + "\t" + overlappingNLRs + "\t" +position + "\t");
				
				if( peptides.containsKey(gene) && peptides.get(gene).containsKey(transcript)) {
					String motifs = peptides.get(gene).get(transcript);
					boolean[] evidence = checkMotifs(motifs);
					
					out.write( evidence[0] + "\t" + evidence[1] + "\t" + evidence[2] + "\t" + motifs);
				
				}else {
					out.write("N/A\tfalse\tN/A\tN/A");
				}
				
				out.newLine();
				
				
				
			}
			
			
			
			
			
		}
		BufferedReader in = new BufferedReader(new FileReader(adr1));

		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			out.write(inputline);
			out.newLine();
		}

		in.close();
		
		
		out.close();
		
		
		
	}
	
	
	
	public static void writeTableS4(
			Hashtable<String, Hashtable<String, int[]>> nlrs,
			Hashtable<String, Hashtable<String, int[]>> genes,
			Hashtable<String, String[]> nlrInfo,
			//Hashtable<String, Hashtable<String, String>> transcripts,
			Hashtable<String, Hashtable<String, String>> peptides,
			File nlrsFasta,
			File manualInspection,
			
			File outputFile)throws IOException{
		
		
		Hashtable<String, String> inspection = new Hashtable<String, String>();
		BufferedReader in = new BufferedReader(new FileReader(manualInspection));
		in.readLine();
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String[] split = inputline.split("\t");
			inspection.put(split[0], split[1]);
			
		}

		in.close();
		
		
		
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		out.write("ID\tOverlapping Genes\tNumber of overlapping genes\tLocus interspesed with \"Ns\"\tOne overlapping gene model considered complete NLR\tLocus position\tStrand\tCompleteness\tMotifs");
		out.newLine();
		
		for(Enumeration<String> myenum1 = nlrs.keys(); myenum1.hasMoreElements();) {
			
			String chromosome = myenum1.nextElement();
			for(Enumeration<String> myenum2 = nlrs.get(chromosome).keys(); myenum2.hasMoreElements();) {
				String nlr = myenum2.nextElement();
				int nlrStart = nlrs.get(chromosome).get(nlr)[0];
				int nlrEnd = nlrs.get(chromosome).get(nlr)[1];
				int nlrStrand = nlrs.get(chromosome).get(nlr)[2];
				
				
				String overlappingGenes = "";
				int numOverlappingGenes = 0;
				boolean locusHasNs = false;
				boolean oneOverlapIsComplete = false;
				String locusPosition = chromosome + ":" + nlrInfo.get(nlr)[3] + "-" + nlrInfo.get(nlr)[4];
				
				
				
				//Check for Ns in the sequence of NLR
				FastaReader fastaReader = new FastaReader(nlrsFasta);
				for (BioSequence seq = fastaReader.readEntry(); seq != null; seq = fastaReader.readEntry()) {
					if( seq.getIdentifier().equalsIgnoreCase(nlr)) {
						if( seq.getSequence().split("N+").length > 1 ) {
							locusHasNs = true;
						}
						
						break;
					}
				}
				fastaReader.close();
				
				
				//Check for overlappingGenes
				HashSet<String> overlappingGeneSet = new HashSet<String>();
				for(Enumeration<String> myenum3 = genes.get(chromosome).keys(); myenum3.hasMoreElements();) {
					String transcript = myenum3.nextElement();
					int transcriptStart = genes.get(chromosome).get(transcript)[0];
					int transcriptEnd = genes.get(chromosome).get(transcript)[1];
					int transcriptStrand = genes.get(chromosome).get(transcript)[2];
					
					if( nlrStrand == transcriptStrand && Math.max(transcriptEnd, nlrEnd) - Math.min(transcriptStart, nlrStart)  < transcriptEnd-transcriptStart + nlrEnd-nlrStart  ) {  //right strand and overlapping
						String gene = transcript.split("\\.")[0];
						overlappingGeneSet.add(gene);
						
						if( peptides.containsKey(gene) && isComplete(peptides.get(gene).get(transcript))) {
							oneOverlapIsComplete = true;
						}
						
					}
				}
				for (Iterator<String> iterator = overlappingGeneSet.iterator(); iterator.hasNext();) {
					String gene = iterator.next();
					overlappingGenes = overlappingGenes + "," + gene;
				}
				numOverlappingGenes = overlappingGeneSet.size();
				if(numOverlappingGenes > 0) {
					overlappingGenes = overlappingGenes.substring(1);
				}
				
				
				
				

				
				out.write(nlr + "\t" + overlappingGenes + "\t" + numOverlappingGenes + "\t" + locusHasNs + "\t" );
				out.write(oneOverlapIsComplete + "\t" + locusPosition + "\t" );
				out.write(nlrInfo.get(nlr)[5] + "\t" + nlrInfo.get(nlr)[2] +"\t" +nlrInfo.get(nlr)[6]);
				
				if( inspection.containsKey(nlr)) {
					out.write("\t"+inspection.get(nlr));
				}
				
				out.newLine();
			}
			
			
			
			
		}
		
		
		
		
		
		out.close();
		
		
	}
	
	private static boolean isComplete(String motifs) {
		boolean[] a = checkMotifs(motifs);
		if(a[0] && a[1] && a[2]) {
			return true;
		}else {
			return false;
		}
		
	}
	
	private static boolean[] checkMotifs(String motifs) {
		
		motifs = "," + motifs + ",";
		
		boolean ploop = false;
		boolean nbarc = false;
		boolean lrr = false;
	
		if(motifs.contains(",1,")){
			ploop = true;
		}	
		if(motifs.contains(",9,") || motifs.contains(",11,") || motifs.contains(",19,") ){
			lrr = true;
		}	
			
		if(motifs.contains(",1,6,4,")
				|| motifs.contains(",1,4,5")
				|| motifs.contains(",6,4,5")
				|| motifs.contains(",4,5,10,")
				|| motifs.contains(",5,10,3,")
				|| motifs.contains(",10,3,12,")
				|| motifs.contains(",3,12,2,")
						){
			nbarc=true;
		}	
		
		
		boolean[] a = {ploop, nbarc, lrr};
		return a;
	}
 	
	
	
	/**
	 * 
	 * read in the NLR-Parser output from CS v1.1 transcripts
	 * @deprecated
	 * @param nlrTransctipts
	 * 		The output File from NLR-Parser ran on CS v1.1 transcripts
	 * @return
	 * 		A Hashtable with gene names as keys. Values are Hashtables, with transcript names as keys and motif sets as values.
	 * @throws IOException
	 */
	public static Hashtable<String, Hashtable<String, String>> readNLRParserTranscripts(File nlrTransctipts)throws IOException{
		//TraesCSU02G059400.1	CNL	complete	forward	42	3354	0	RSIVRQICVNSLQETGEEGKTTIGAQVLKMMSMTKEDGLACEFKRYLNDKSYLIVLNGIHTTEEWDCIKTCFPTNKKGSRIIVCTKQVEVACLCIEAEDEALVHKMLFADQSLYAFYSKERRNSKEKGSIPFVVSTGVNSSTHKNILTQTEAMATLEESRLIGRGNEKEEIIKLISNKDLRQFHVISLWGMGGIGKTTLVKDIYQSQENSS	false	17,16,6,1,4,5,1,6,4,5,10,3,12,2,8,7,11,9,11
		Hashtable<String, Hashtable<String, String>> h = new Hashtable<String, Hashtable<String, String>>();
		
		BufferedReader in = new BufferedReader(new FileReader(nlrTransctipts));
		in.readLine();
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {

			String[] split = inputline.split("\t");
			
			if( !split[3].equalsIgnoreCase("forward")) { //if the nlr is called on the reverse strand of the gene the gene has been annotated wrongly.
				continue;
			}
			
			String gene = split[0].split("\\.")[0];
			
			if(!h.containsKey(gene)) {
				h.put(gene, new Hashtable<String, String>() );
			}
			
			h.get(gene).put(split[0],split[9]);
			
		}

		in.close();
		return h;
		
	}
	
	/**
	 * 
	 * read in the NLR-Parser output from CS v1.1 peptides
	 * This works with NLR-Parser v. 1.0
	 * 
	 * @param nlrTransctipts
	 * 		The output File from NLR-Parser ran on CS v1.1 peptides
	 * @return
	 * 		A Hashtable with gene names as keys. Values are Hashtables, with transcript names as keys and motif sets as values.
	 * @throws IOException
	 */
	public static Hashtable<String, Hashtable<String, String>> readNLRParserPeptides(File nlrPeptides)throws IOException{
		//TraesCSU02G059400.1	CNL	complete	forward	42	3354	0	RSIVRQICVNSLQETGEEGKTTIGAQVLKMMSMTKEDGLACEFKRYLNDKSYLIVLNGIHTTEEWDCIKTCFPTNKKGSRIIVCTKQVEVACLCIEAEDEALVHKMLFADQSLYAFYSKERRNSKEKGSIPFVVSTGVNSSTHKNILTQTEAMATLEESRLIGRGNEKEEIIKLISNKDLRQFHVISLWGMGGIGKTTLVKDIYQSQENSS	false	17,16,6,1,4,5,1,6,4,5,10,3,12,2,8,7,11,9,11
		//TraesCS6A02G290500.1	N/A	partial	TraesCS6A02G290500.1:190	TraesCS6A02G290500.1:427	N/A	false	11,12,11,12,9,11
		Hashtable<String, Hashtable<String, String>> h = new Hashtable<String, Hashtable<String, String>>();
		
		BufferedReader in = new BufferedReader(new FileReader(nlrPeptides));
		in.readLine();
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {

			String[] split = inputline.split("\t");
			
			
			String gene = split[0].split("\\.")[0];
			
			if(!h.containsKey(gene)) {
				h.put(gene, new Hashtable<String, String>() );
			}
			
			h.get(gene).put(split[0],split[7]);
			
		}

		in.close();
		return h;
		
	}
	
	
	/**
	 * 
	 * Read in the output of NLR-Annotator (tsv format) and converts this to a Hashtable.
	 * 
	 * @param nlrTxt
	 * 			NLR-Annotator output in tsv format
	 * @return
	 * 			A Hashtable with NLR loci as keys and the entire tsv entry as value
	 * @throws IOException
	 */
	public static Hashtable<String, String[]> getNLRInfo(File nlrTxt)throws IOException{
		Hashtable<String, String[]> h = new Hashtable<String, String[]>();
		
		BufferedReader in = new BufferedReader(new FileReader(nlrTxt));

		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String[] split = inputline.split("\t");
			h.put(split[1], split);
		}

		in.close();
		
		return h;
	}
	
	
	
	
	public static Hashtable<String, Hashtable<String, int[]>> getNLRs(File file)throws IOException{
		Hashtable<String, Hashtable<String, int[]>> h = new Hashtable<String, Hashtable<String, int[]>>();
		
		
		BufferedReader in = new BufferedReader(new FileReader(file));

		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			//chr2D	chr2D_nlr_1	complete	1602055	1609183	+	17,16,1,6,4,3,2,11,9,9
			String[] split = inputline.split("\t");
			//if(!split[2].equalsIgnoreCase("complete")) {
				//continue;
			//}
			String id = split[1];
			int start = Integer.parseInt(split[3]);
			int end = Integer.parseInt(split[4]);
			int strand = 1;
			if(split[5].equalsIgnoreCase("-")) {
				strand = -1;
			}
			
			if( !h.containsKey(split[0])) {
				h.put(split[0], new Hashtable<String, int[]>());
			}
			int[] a = {start, end, strand};
			h.get(split[0]).put(id, a);

		}

		in.close();
		
		return h;
	}
	
		
	
	public static Hashtable<String, Hashtable<String, int[]>> getGFF(File hcGFF, File lcGFF)throws IOException{
		
		Hashtable<String, Hashtable<String, int[]>> h = new Hashtable<String, Hashtable<String, int[]>>();
		
		File[] gffs = {hcGFF, lcGFF};
		for( int index = 0; index < gffs.length; index ++) {
			
			BufferedReader in = new BufferedReader(new FileReader(gffs[index]));

			for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
				//chr1A	IWGSC_March2017	gene	41202	41522	31	+	.	ID=TraesCS1A01G000100LC;primconf=LC
				
				if(inputline.startsWith("#")) {
					continue;
				}
				
				String[] split  = inputline.split("\t");
				
				if(!split[2].equalsIgnoreCase("mRNA")) {
					continue;
				}
				
				String id = split[8].split("ID=")[1].split(";")[0];
				int start = Integer.parseInt(split[3]);
				int end = Integer.parseInt(split[4]);
				int strand = 1;
				if(split[6].equalsIgnoreCase("-")) {
					strand = -1;
				}
				
				if( !h.containsKey(split[0])) {
					h.put(split[0], new Hashtable<String, int[]>());
				}
				int[] a = {start, end, strand};
				h.get(split[0]).put(id, a);
				
			}

			in.close();
			
		}
		
		return h;
		
	}
	
	
}
