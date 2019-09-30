

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Vector;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import support.BioSequence;
import support.FastaReader;

public class PreparePhylogenetics {

	public static void main(String[] args) {
		try {
			
			/*
			 * ******** Input Data ************
			 */
			
			// IWGSC gene annotation v. 1.1
			File hcCDS = new File("steuernb/nlr_paper/v1.1_analysis/inputData/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_cds.fasta");
			File lcCDS = new File("steuernb/nlr_paper/v1.1_analysis/inputData/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_LC_20170706_cds.fasta");
			File hcPEP = new File("steuernb/nlr_paper/v1.1_analysis/inputData/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_pep.fasta");
			File lcPEP = new File("steuernb/nlr_paper/v1.1_analysis/inputData/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_LC_20170706_pep.fasta");
			File hcGFF = new File("steuernb/nlr_paper/v1.1_analysis/inputData/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.gff3");
			File lcGFF = new File("steuernb/nlr_paper/v1.1_analysis/inputData/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_LC_20170706.gff3");
			
			
			File rgenesPEP =  new File("steuernb/nlr_paper/v1.1_analysis/inputData/Rgenes.protein.fasta");
			
			
			
			// HMMScan
			File hmm = new File("steuernb/nlr_paper/v1.1_analysis/HMMScan/CS_1.1_pep.hmmscan.txt");
			File hmmRgene = new File("steuernb/nlr_paper/v1.1_analysis/HMMScan/Rgenes.protein.hmmscan.txt");
			
			
			// MAST on NLR peptide sequences
			File mastPepFile = new File("steuernb/nlr_paper/v1.1_analysis/TaeCS_NlrPep.mast.xml");
			File mastPepRgeneFile = new File("steuernb/nlr_paper/v1.1_analysis/Rgenes_mast.xml");
			
			// transcripts with complete NLRs selected from TableS5 This is also Table S6
			File nlrTranscriptsFile = new File("steuernb/nlr_paper/v1.1_analysis/intermediate/NLR-transcripts.txt");
			
			// transcripts from the ID clade
			File nlrIDClade = new File("steuernb/nlr_paper/v1.1_analysis/multiple_alignment_NB-ARC/ID_clade_member.txt");
			File connectorClade = new File("steuernb/nlr_paper/v1.1_analysis/multiple_alignment_NB-ARC/Connector_clade_member.txt");
			
			
			
			/*
			 * ******** Output Data ************
			 */
			File nbarc_sequenceFile = new File("steuernb/nlr_paper/v1.1_analysis/TaeCS_NB-ARC.fasta");
			File pep_sequenceFile = new File("steuernb/nlr_paper/v1.1_analysis/TaeCS_NlrPep.fasta");
			File iTOL_domainsFile =  new File("steuernb/nlr_paper/v1.1_analysis/multiple_alignment_NB-ARC/iTOL_domanis.txt");
			File iTOL_connectors =  new File("steuernb/nlr_paper/v1.1_analysis/multiple_alignment_NB-ARC/iTOL_connectors.txt");
			File iTOL_introns =  new File("steuernb/nlr_paper/v1.1_analysis/multiple_alignment_NB-ARC/iTOL_introns.txt");
			File tableS7 = new File("steuernb/nlr_paper/v1.1_analysis/outputData/tableS7.txt");
			File tableS6 = new File("steuernb/nlr_paper/v1.1_analysis/outputData/tableS7.txt");
			
			/* 
			 ************************************** 
			 */
				
			
			
			
			HashSet<String> nlrTranscripts = readNLRTranscriptList(nlrTranscriptsFile, hcCDS, lcCDS);
			//System.out.println(nlrTranscripts.size());
			
			/*
			 * extract the NB-ARC domains from aa sequences of NLR genes
			 */
			 
			extractNBARC(nlrTranscripts, hmm,hmmRgene, hcPEP, lcPEP, rgenesPEP, nbarc_sequenceFile);	
			
			
			
			/*
			 * just write the aa sequences for all nlr genes. This is input for MAST ( source meme-4.9.1; mast ~/meme.xml TaeCS_NlrPep.fasta ) 
			 */
			
			extractNlrPepSequences(nlrTranscripts, hcPEP, lcPEP,rgenesPEP, pep_sequenceFile);	
			
			
			/*
			 * Generate the motif and domain annotation for iTOL 
			 */
			
			HashMap<String, Vector<String>> domains =  getDomains(nlrTranscripts,  hmm, hmmRgene, 1E-2);
			HashMap<String, Vector<String>> motifs  =  getMotifs (nlrTranscripts,  mastPepFile, mastPepRgeneFile );
			
			writeITOLDomains( domains,  motifs, iTOL_domainsFile,  tableS7);
			
			
			
			/*
			 * Write connectors
			 */
			
			writeITOLTandemPairs(nlrTranscripts, hcGFF, lcGFF, 20000, nlrIDClade, connectorClade, iTOL_connectors);
			
			
			/*
			 * Write Intron/Exon structure along with motifs
			 */
			writeIntronExonForITOL2(nlrTranscripts, motifs, hcGFF, lcGFF, iTOL_introns);
			
		} catch (Exception e) {
			e.printStackTrace();
		}

	}
	
	public static void writeTableS6(HashMap<String, Vector<String>> domains, File tableS6)throws IOException{
		
		for(Iterator<String> iterator = domains.keySet().iterator(); iterator.hasNext();) {
			String gene = iterator.next();
			
			System.out.println(gene);
			
		}
		
	}
	
	
	public static HashSet<String> readNLRTranscriptList(File nlrTranscriptsFile)throws IOException{
		System.out.println("reading NLR Transcript list " + nlrTranscriptsFile.getAbsolutePath());
		HashSet<String> h = new HashSet<String>();
		BufferedReader in = new BufferedReader(new FileReader(nlrTranscriptsFile));

		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			h.add(inputline.split("\t")[0].trim());
		}

		in.close();
		return h;
	}
	
	/**
	 * 
	 * Get one transcript per gene. Always the longest.
	 * 
	 * @param nlrTranscriptsFile
	 * @param inputCDS
	 * @return
	 * @throws IOException
	 */
	public static HashSet<String> readNLRTranscriptList(File nlrTranscriptsFile, File hcCDS, File lcCDS)throws IOException{
		System.out.println("reading NLR Transcript list " + nlrTranscriptsFile.getAbsolutePath());
		
		
		HashSet<String> h = new HashSet<String>();
		
		BufferedReader in = new BufferedReader(new FileReader(nlrTranscriptsFile));

		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			h.add(inputline.split("\t")[0].trim());
		}

		in.close();
		
		
		HashMap<String, String> geneToTranscript = new HashMap<String, String>();
		HashMap<String, Integer> transcriptLength = new HashMap<String, Integer>();
		
		File[] cds = {hcCDS, lcCDS};
		for( int i = 0; i< cds.length; i++) {
			FastaReader fastaReader = new FastaReader(cds[i]);
			for (BioSequence seq = fastaReader.readEntry(); seq != null; seq = fastaReader.readEntry()) {
				String gene = seq.getIdentifier().split("\\.")[0];
				
				if( h.contains(seq.getIdentifier())) {
					
					if( !geneToTranscript.containsKey(gene) || transcriptLength.get(gene) < seq.getLength()  ) {
						geneToTranscript.put(gene, seq.getIdentifier());
						transcriptLength.put(gene, seq.getLength());
					}
					
				}
				
			}
			fastaReader.close();
		}
		
		
		
		h.clear();
		h.addAll(geneToTranscript.values());
		
		return h;
	}
	
	public static void extractNlrPepSequences(HashSet<String> nlrTranscripts, File hcPEP, File lcPEP,File rgenesPEP, File outputFile)throws IOException{
		
		
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		
		
		File[] peps = {hcPEP, lcPEP, rgenesPEP};
		for( int i = 0; i< peps.length; i++) {
			FastaReader fastaReader = new FastaReader(peps[i]);
			for (BioSequence seq = fastaReader.readEntry(); seq != null; seq = fastaReader.readEntry()) {
				if( nlrTranscripts.contains(seq.getIdentifier()) || i==2) {
					out.write(seq.getFastaString(100));
				}
			}
			fastaReader.close();
		}
		out.close();
	}
	
	
	
	/**
	 * 
	 * Extract the NB-ARC domain sequence for running phylogenetics
	 * 
	 * @param nlrTranscripts
	 * 				list of IDs that are included
	 * @param hmm
	 * 				Output file from hmmscan of wheat protein sequences
	 * @param hmmRgene
	 * 				Output file from hmmscan on protein sequences from R genes that are to be included in the phylogenetics
	 * @param hcPEP
	 * 				Protein sequences from high confidence genes from IWGSC 1.1 gene annotation on Chinese Spring
	 * @param lcPEP
	 * 				Protein sequences from low confidence genes from IWGSC 1.1 gene annotation on Chinese Spring
	 * @param rgenePEP
	 * 				Protein sequences from R genes
	 * @param outputFile
	 * 				location of output file
	 * @throws IOException
	 */
	public static void extractNBARC(HashSet<String> nlrTranscripts, File hmm,File hmmRgene, File hcPEP, File lcPEP, File rgenePEP, File outputFile)throws IOException{
		
		System.out.println("extracting NB-ARC domains");
		
		HashMap<String, BioSequence> genes = new HashMap<String, BioSequence>();
		File[] peps = { hcPEP, lcPEP, rgenePEP };
		for( int i = 0; i< peps.length; i++) {
			FastaReader fastaReader = new FastaReader(peps[i]);
			for (BioSequence seq = fastaReader.readEntry(); seq != null; seq = fastaReader.readEntry()) {
				if( nlrTranscripts.contains(seq.getIdentifier()) || i == 2) {
					genes.put(seq.getIdentifier(), seq);
				}
			}
			fastaReader.close();
		}
		
		
		
		
		
		HashMap<String, int[]> nbarcs = new HashMap<String, int[]>();
		
		File[] hmms = {hmm, hmmRgene};
		for( int i = 0; i< hmms.length; i++) {
			BufferedReader in = new BufferedReader(new FileReader(hmms[i]));

			for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
				//NB-ARC               PF00931.21   288 TraesCS1A02G003546.1 -            368   1.6e-55  188.1   0.1   1   1     1e-58   2.1e-55  187.7   0.1    16   287    92   364    80   365 0.89 NB-ARC domain
				String[] split = inputline.split("\\s+");
				if( split[0].equalsIgnoreCase("NB-ARC") && genes.containsKey(split[3]) ) {
					
					String gene = split[3];
					
					int start = Integer.parseInt(split[19]);
					int end = Integer.parseInt(split[20]);
					
					int[] a = {0,0};
					if( nbarcs.containsKey(gene)) { //there might be more than one NB-ARC in the gene. We take the longest for now.
						a = nbarcs.get(gene);
					}
					if( end-start > a[1] - a[0] ) {
						a[0] = start;
						a[1] = end;
					}
					
					
					nbarcs.put(gene, a);
					
				}
				
			}
			in.close();	
		}
		
		

		
			
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		for(Iterator<String> iterator = genes.keySet().iterator(); iterator.hasNext();) {
			String gene = iterator.next();
			
			if( nbarcs.containsKey(gene)) {
				int[] a = nbarcs.get(gene);
				
				//String nbarc = genes.get(gene).getSequence().substring( (a[0]-1)*3 , a[1] * 3); //convert aa pos to nucleotide pos   //deprecated. I need 
				String nbarc = genes.get(gene).getSequence().substring( a[0] , a[1] ); 
				BioSequence seq = new BioSequence(gene, nbarc);
				seq.setDescription( "NBARC, " + a[0] + "-" + a[1] );
				out.write(seq.getFastaString());
				
			}else {
				System.out.println("no NB-ARC domain found for " +gene);
			}
			
			
		}
		
		out.close();
		
		
		
	}
	
	
	
	/**
	 * 
	 * Extract additional integrated domains from hmmscan
	 * 
	 * @param nlrTranscriptNames
	 * 			list of IDs that are included in analysis
	 * @param hmm
	 * 			output of hmmscan
	 * @param hmmRgene
	 * 			output of hmmscan on R genes
	 * @param evalueThreshold
	 * 			evalue threshols for hmmscan output
	 * @return
	 * @throws IOException
	 */
	public static HashMap<String, Vector<String>> getDomains(HashSet<String> nlrTranscriptNames, File hmm,File hmmRgene, double evalueThreshold)throws IOException{
		
		
		
		HashMap<String, Vector<String> > h = new HashMap<String, Vector<String> >();
		
		for(Iterator<String> iterator = nlrTranscriptNames.iterator(); iterator.hasNext();){
			h.put(iterator.next(), new Vector<String> ());
		}
		
		File[] files = {hmm, hmmRgene};
		
		for( int i = 0; i<files.length; i++) {
			BufferedReader in = new BufferedReader(new FileReader(files[i]));

			for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
				if( inputline.startsWith("#")){
					continue;
				}
				String[] split = inputline.split("\\s+");
				String id = split[3];
				
				if(i==1 && !h.containsKey(id)) {
					h.put(id, new Vector<String>() );
				}
				
				if(h.containsKey(id) ){
					double evalue = Double.parseDouble(split[6]);
					if(evalue <= evalueThreshold){
						h.get(id).add(inputline);
					}
				}
			}

			in.close();
		}
		
	
		
		for(Iterator<String> iterator = h.keySet().iterator(); iterator.hasNext();){
			String id = iterator.next();
			Vector<String> v = h.get(id);
			Collections.sort(v, new Comparator<String>(){public int compare(String s1, String s2){
															String[] split1 = s1.split("\\s+");
															String[] split2 = s2.split("\\s+");
															if(Double.parseDouble(split1[6]) < Double.parseDouble(split2[6])){
																return -1;
															}
															if(Double.parseDouble(split1[6]) > Double.parseDouble(split2[6])){
																return 1;
															}
															return 0;}} );
			Vector<String> v2 = new Vector<String>();
			for(Iterator<String> iterator2 = v.iterator(); iterator2.hasNext();){
				String s1 = iterator2.next();
				boolean good = true;
				for(Iterator<String>iterator3 = v2.iterator(); iterator3.hasNext();){
					String s2 = iterator3.next();
					String[] split1 = s1.split("\\s+");
					String[] split2 = s2.split("\\s+");
					int start1 = Integer.parseInt(split1[19]);
					int end1 = Integer.parseInt(split1[20]);
					int start2 = Integer.parseInt(split2[19]);
					int end2 = Integer.parseInt(split2[20]);
					if(end1-start1 + end2-start2 >= Math.max(end1, end2) -Math.min(start1, start2)){
						good = false;
						break;
					}
				}
				if(good){
					v2.add(s1);
				}
			}
			
			h.put(id, v2);
			
		}
		
		
		
		
		return h;
	}
	
	
	/**
	 * 
	 * Extract the position of motifs found by MAST (using the NLR-Parser meme.xml)
	 * 
	 * 
	 * @param nlrNames
	 * 			List of IDs to be included in analysis
	 * @param mastPepFile
	 * 			xml output from MAST on IDs
	 * @param mastPepRgeneFile
	 * 			xml output from MAST on Rgenes
	 * @return
	 * 		HashMap has IDs as keys. The Values are Strings that are components for the iTOL files.
	 * 
	 * @throws ParserConfigurationException
	 * @throws IOException
	 * @throws SAXException
	 */
	public static HashMap<String, Vector<String>> getMotifs(HashSet<String> nlrNames, File mastPepFile, File mastPepRgeneFile )throws ParserConfigurationException , IOException, SAXException{
		
		HashMap<String,Vector<String>> h = new HashMap<String,Vector<String>>();
		for(Iterator<String> iterator = nlrNames.iterator(); iterator.hasNext();){
			h.put(iterator.next(), new Vector<String> ());
		}
		
		File[] files = {mastPepFile, mastPepRgeneFile};
		
		for( int i = 0; i< files.length; i++) {
			DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
			DocumentBuilder db = dbf.newDocumentBuilder();
			Document dom = db.parse(files[i]);
			
				
			Element rootElement = dom.getDocumentElement();
			NodeList nodeList1 = rootElement.getElementsByTagName("sequences");
			for( int index1 = 0; index1< nodeList1.getLength(); index1++){
				Element sequences = (Element) nodeList1.item(index1);
				
				NodeList nodeList2 = sequences.getElementsByTagName("sequence");
				for( int index2 =0; index2 < nodeList2.getLength(); index2++){
					Element sequence = (Element) nodeList2.item(index2);
					
					String id = sequence.getAttribute("name");
					
					if( i==1 && !h.containsKey(id)) {
						h.put(id, new Vector<String>());
					}
					
					if(!h.containsKey(id)){
						continue;
					}
					
					
					NodeList nodeList3 = sequence.getElementsByTagName("hit");
					for( int index3 = 0; index3< nodeList3.getLength(); index3++){
						Element hit = (Element)nodeList3.item(index3);
						int pos = Integer.parseInt(hit.getAttribute("pos"));
						int end = pos + hit.getAttribute("match").length();
						String name = hit.getAttribute("motif");
						String color = getMotifColor(name);
						double pvalue = Double.parseDouble(hit.getAttribute("pvalue"));
						
						if(color.equalsIgnoreCase("none") || pvalue > 1E-5){
							continue;
						}
						
						h.get(id).add("RE|"+pos+"|"+end+"|"+color+"|"+name);
						
					}
					
				}
				
			}
		}
			
	
		return h;	
	}		
	
	public static void writeITOLDomains(HashMap<String,Vector<String>> domains, HashMap<String,Vector<String>> motifs, File outputFile, File tableS7)throws IOException{
		
		HashMap<String, HashSet<String>> domainCount = new HashMap<String, HashSet<String>>();
		HashMap<String, HashSet<String>> domainCount5 = new HashMap<String, HashSet<String>>();
		HashMap<String, HashSet<String>> domainCount10 = new HashMap<String, HashSet<String>>();
		
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		out.write("DATASET_DOMAINS"); out.newLine();
		out.write("SEPARATOR COMMA"); out.newLine();
		out.write("DATASET_LABEL,Domains"); out.newLine();
		out.write("COLOR,#ff0000"); out.newLine();
		out.write("DATASET_SCALE,100-amino_acid_100-#0000ff,500-line at aa500-#ff0000,1000-aa1000-#00ff00"); out.newLine();
		out.write("DATA"); out.newLine();
		
		
		
		
		
		for(Iterator<String> iterator = motifs.keySet().iterator(); iterator.hasNext();){
			String id = iterator.next();
			Vector<String> motifList = motifs.get(id);
			Vector<String> domainList = domains.get(id);
			int length = Integer.parseInt(domainList.iterator().next().split("\\s+")[5]);
			out.write(id + ","+length);
			
			for(Iterator<String> iterator2 = domainList.iterator(); iterator2.hasNext();){
				String[] split = iterator2.next().split("\\s+");
				//t(",RE|"+pos+"|"+end+"|"+color+"|"+name);
				
				if(split[0].startsWith("NB-ARC") || split[0].startsWith("LRR") || split[0].startsWith("AAA")|| split[0].startsWith("NACHT")){
					continue;
				}
				
				
				
				double evalue = Double.parseDouble(split[6]);
				
				
				out.write(",RE|"+Integer.parseInt(split[19]) + "|" +Integer.parseInt(split[20]) + "|" + getDomainColor(split[0]) + "|"+split[0]);
				
				if(!domainCount.containsKey(split[0])) {
					domainCount.put(split[0], new HashSet<String>());
					domainCount10.put(split[0], new HashSet<String>());
					domainCount5.put(split[0], new HashSet<String>());
				}
				
				domainCount.get(split[0]).add(id);
				if(evalue <= 1E-5) {domainCount5.get(split[0]).add(id);}
				if(evalue <= 1E-10) {domainCount10.get(split[0]).add(id);}
			}
			
			for(Iterator<String> iterator2 = motifList.iterator(); iterator2.hasNext();){
				out.write(","+iterator2.next());
			}
			
			out.newLine();
		}
		
		out.close();
		
		 out = new BufferedWriter(new FileWriter(tableS7));
		
		 HashSet<String>transcriptCount = new HashSet<String>();
		 HashSet<String>transcriptCount5 = new HashSet<String>();
		 HashSet<String>transcriptCount10 = new HashSet<String>();
		 
		 
		 for(Iterator<String> iterator = domainCount.keySet().iterator(); iterator.hasNext();) {
			 String key = iterator.next();
			 int d1 = domainCount.get(key).size();
			 int d5 = domainCount5.get(key).size();
			 int d10 = domainCount10.get(key).size();
			 
			 out.write(key + "\t" + d1 + "\t" + d5 + "\t" + d10);
			 
			 
			 transcriptCount.addAll(domainCount.get(key));
			 transcriptCount5.addAll(domainCount5.get(key));
			 transcriptCount10.addAll(domainCount10.get(key));
			 
			 out.newLine();
		 }
		out.close();
		
		System.out.println("Number of genes with ID: " + transcriptCount.size() + "\t" + transcriptCount5.size() + "\t" + transcriptCount10.size());
	}
		
	
	
	
	
	private static String getDomainColor(String domain){
		
		if(domain.equalsIgnoreCase("Pkinase")){return "#ff0010";}
		if(domain.equalsIgnoreCase("DDE_Tnp_4")){return "#800080";}
		if(domain.equalsIgnoreCase("zf-BED")){return "#f0a3ff";}
		if(domain.startsWith("Kelch")){return "#ffa8bb";}
		if(domain.equalsIgnoreCase("Motile_Sperm")){return "#c20088";}
		if(domain.equalsIgnoreCase("Jacalin")){return "#990000";}
		if(domain.equalsIgnoreCase("Pkinase_Tyr")){return "#de0000";}
		
		
		return "#8f7c00";
	}
	
	private static String getMotifColor(String name){
		if(name.equalsIgnoreCase("motif_1")){return "#00ffff";}
		if(name.equalsIgnoreCase("motif_6")){return "#00e1ff";}
		if(name.equalsIgnoreCase("motif_4")){return "#00c3ff";}
		if(name.equalsIgnoreCase("motif_5")){return "#00a5ff";}
		if(name.equalsIgnoreCase("motif_10")){return "#0087ff";}
		if(name.equalsIgnoreCase("motif_3")){return "#0069ff";}
		if(name.equalsIgnoreCase("motif_12")){return "#004bff";}
		if(name.equalsIgnoreCase("motif_2")){return "#002dff";}
		if(name.equalsIgnoreCase("motif_9")){return "#008000";}
		if(name.equalsIgnoreCase("motif_11")){return "#00c000";}
		if(name.equalsIgnoreCase("motif_19")){return "#00ff00";}
		if(name.equalsIgnoreCase("motif_17")){return "#993f00";}
		if(name.equalsIgnoreCase("motif_16")){return "#ffa405";}
	
		return "none";
	
	}
	
	
	public static void writeITOLTandemPairs(HashSet<String> nlrTranscripts, File hcGFF, File lcGFF, int maxDistance, File nlrIDClade,File connectorClade, File outputFile )throws IOException{
		
		/*
		 *read in all genes from gffs.  
		 */
		File[] gffs = {hcGFF, lcGFF};
		Hashtable<String, int[]> genePositions = new Hashtable<String, int[]>();
		
		Vector<String> geneOrder = new Vector<String>(); //this is to check how many genes are between an NLR pair.
				
		
		for( int i = 0; i< gffs.length; i++) {
			File gff = gffs[i];
			
			BufferedReader in = new BufferedReader(new FileReader(gff));

			for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
				if(inputline.startsWith("#")) {
					continue;
				}
				
				String[] split = inputline.split("\t");
				
				if( split[2].equalsIgnoreCase("gene")) {
					String id = split[8].split("ID=")[1].split(";")[0];
					geneOrder.add(id + "\t" + split[0] + "\t" + split[3]); //add each gene to this vector together with chr and position
				}
				
				if(!split[2].equalsIgnoreCase("mRNA")){
					continue;
				}
				String id = inputline.split("ID=")[1].split(";")[0];
				if(!nlrTranscripts.contains(id)){
					continue;
				}
				int direction = 0;
				if(split[6].equalsIgnoreCase("+")){
					direction = 1;
				}
				
				int[] a = {Integer.parseInt(split[3]), Integer.parseInt(split[4]),direction}; //record start, end and direction for each NLR 
				genePositions.put(id, a);
				
				
			}

			in.close();
		}
		
		
		//order the geneOrder
		Collections.sort(geneOrder, new Comparator<String>(){ public int compare(String s1, String s2) {
																			String[] split1 = s1.split("\t");
																			String[] split2 = s2.split("\t");
																			int myreturn = 0;
																			if(split1[1].equalsIgnoreCase(split2[1])) {
																				int i1 = Integer.parseInt(split1[2]);
																				int i2 = Integer.parseInt(split2[2]);
																				if( i1< i2) {myreturn =  -1;}
																				if( i1> i2) {myreturn =   1;}
																				if( i1==i2) {myreturn =  0;}
																				
																				//System.out.println(myreturn);
																			}else {
																				myreturn =  split1[1].compareTo(split2[1]);
																			}
																			
																			return myreturn;
																			}});
		
		
		
		
		// load the IDs from NLR clade
		HashSet<String> idClade = new HashSet<String>();
		BufferedReader in = new BufferedReader(new FileReader(nlrIDClade));
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			idClade.add(inputline.trim());
		}
		in.close();
		
		//load the IDs from Connector clade
		HashSet<String> conClade = new HashSet<String>();
		 in = new BufferedReader(new FileReader(connectorClade));
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			conClade.add(inputline.trim());
		}
		in.close();
		
		
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		out.write("DATASET_CONNECTION"); out.newLine();
		out.write("SEPARATOR COMMA"); out.newLine();
		out.write("DATASET_LABEL,NLR_pairs"); out.newLine();
		out.write("COLOR,#ff0ff0"); out.newLine();
		out.write("DRAW_ARROWS,0"); out.newLine();
		out.write("ARROW_SIZE,20"); out.newLine();
		out.write("MAXIMUM_LINE_WIDTH,1"); out.newLine();
		out.write("CURVE_ANGLE,0"); out.newLine();
		out.write("CENTER_CURVES,1"); out.newLine();
		out.write("ALIGN_TO_LABELS,1"); out.newLine();
		out.write("DATA"); out.newLine();
		
		Hashtable<String, int[]> connections = new Hashtable<String, int[]> ();
		HashSet<String> wasThere = new HashSet<String>();
		/*
		 * now go through each combination of genes and check their direction and distance
		 */
		for(Enumeration<String> myenum1 = genePositions.keys(); myenum1.hasMoreElements();) {
			String key = myenum1.nextElement();
			
			for(Enumeration<String> myenum2 =genePositions.keys(); myenum2.hasMoreElements();) {
				String partner = myenum2.nextElement();
				
				if(key.equalsIgnoreCase(partner)) {
					continue;
				}
				
				//TraesCS1A
				if( ! key.substring(0, 9).equalsIgnoreCase(partner.substring(0, 9))) { //check if genes are on the same chromosome. chromosome is encoded in gene name.
					continue;
				}
				
				int[] a = genePositions.get(key);
				int[] b = genePositions.get(partner);
				
				if( a[2] != b[2] ) {
					int distance = Math.max(a[0], b[1]) - Math.min(a[0], b[1]);
					if( distance <= maxDistance) {
						
						
						
						String color = "#EBEBEB";
						if(distance<20000) { color = "#B7FFBB"; }                  // "#929292"; }
						if(distance< 5000) { color = "#339A3A"; }                   //"#0D0D0D"; }
						
						
						String[] v = {key, partner};
						if( b[0] < a[0]) {  //have a pair of NLRs in alphabetical order to not report them twice
							v[0] = partner;
							v[1] = key;
						}
						
						if(wasThere.contains(v[0] + " "+ v[1])) {
							continue;
						}
						wasThere.add(v[0] + " " + v[1]);
						
						String s = v[0] + "," + v[1] + ",1," + color+",normal";
						
						out.write(s);
						out.newLine();
						
						int[] c = {distance,0,0, 0};
						if( idClade.contains(v[0]) ) {c[1] = 1;} 
						if( idClade.contains(v[1]) ) {c[1] = 2;}
						if(conClade.contains(v[0]) ) {c[2] = 1;}
						if(conClade.contains(v[1]) ) {c[2] = 2;}	
							
						String id1 = 	v[0].split("\\.")[0];
						String id2 =    v[1].split("\\.")[0];
						
						boolean found = false;
						int numGenes = 0;
						for(Enumeration<String> myenum = geneOrder.elements(); myenum.hasMoreElements();) {
							String[] split = myenum.nextElement().split("\t");
							if(split[0].equalsIgnoreCase(id2)) {
								break;
							}
							if(found) {
								numGenes++;
							}
							if(split[0].equalsIgnoreCase(id1)) {
								found = true;
							}
							
						}
						c[3] = numGenes;
						
						connections.put( id1 + "\t" +id2  , c);
						
						
					}
					  
				}
				
				
			}
			
			
		}
		
		out.close();
		
		
		for(Enumeration<String> myenum = connections.keys(); myenum.hasMoreElements();) {
			String key = myenum.nextElement();
			
			System.out.println( key+ "\t"+connections.get(key)[0] + "\t" + connections.get(key)[1]+ "\t" + connections.get(key)[2] + "\t" + connections.get(key)[3]);
			
		}
		
		
		
		
	}
	
	
	
	public static void writeIntronExonForITOL(HashSet<String> transcripts, HashMap<String,Vector<String>> motifs, File gffHC, File gffLC, File outputFile )throws IOException{
		
		Hashtable<String, Vector<int[]>> exons = new Hashtable<String, Vector<int[]>>();
		Hashtable<String, int[]> utr3 = new Hashtable<String,int[]>();
		Hashtable<String, int[]> utr5 = new Hashtable<String,int[]>();
		
		File[] gffs = {gffHC, gffLC};
		for( int i = 0; i< gffs.length; i++) {
			File gff = gffs[i];
			BufferedReader in = new BufferedReader(new FileReader(gff));

			for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
				if(inputline.startsWith("#")) {
					continue;
				}
				String[] split = inputline.split("\t");
				//chr1A	IWGSC_v1.1_201706	exon	40098	40731	.	-	.	ID=TraesCS1A02G000100.1.exon1;Parent=TraesCS1A02G000100.1
				if( split[2].equalsIgnoreCase("exon")) {
					String transcript = split[8].split("Parent=")[1].split(";")[0];
					if(transcripts.contains(transcript)) {
						int start = Integer.parseInt(split[3]);
						int end = Integer.parseInt(split[4]);
						int strand = 1;if(split[6].equalsIgnoreCase("-")) {	strand  = -1;	}
						
						if( !exons.containsKey(transcript)) {
							exons.put(transcript, new Vector<int[]>());
						}
						int[] a = {start, end, strand};
						exons.get(transcript).addElement(a);
					}
					
				}
				
				if( split[2].equalsIgnoreCase("three_prime_UTR")) {
					String transcript = split[8].split("Parent=")[1].split(";")[0];
					if(transcripts.contains(transcript)) {
						int start = Integer.parseInt(split[3]);
						int end = Integer.parseInt(split[4]);
						int strand = 1;
						if(split[6].equalsIgnoreCase("-")) {
							strand  = -1;
						}
						int[] a = {start,end, strand};
						utr3.put(transcript, a);
					}	
				}
				if( split[2].equalsIgnoreCase("five_prime_UTR")) {
					String transcript = split[8].split("Parent=")[1].split(";")[0];
					if(transcripts.contains(transcript)) {
						int start = Integer.parseInt(split[3]);
						int end = Integer.parseInt(split[4]);
						int strand = 1;
						if(split[6].equalsIgnoreCase("-")) {
							strand  = -1;
						}
						int[] a = {start,end, strand};
						utr5.put(transcript, a);
					}	
				}
				
			}

			in.close();
			
		}
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		out.write("DATASET_DOMAINS");out.newLine();
		out.write("SEPARATOR COMMA");out.newLine();
		out.write("DATASET_LABEL,Introns + Domains");out.newLine();
		out.write("COLOR,#00ff00");out.newLine();
		out.write("DATASET_SCALE,1000-1000-#0000ff,2000-2000-#ff0000,3000-1000-#00ff00");out.newLine();
		out.write("DATA");out.newLine();
		
		
		
		for( Enumeration<String> myenum = exons.keys(); myenum.hasMoreElements();) {
			String transcript = myenum.nextElement();
			Vector<int[]> exonVector = exons.get(transcript);
			Collections.sort(exonVector, new Comparator<int[]>() {public int compare(int[] a1, int[] a2) {
																if( a1[0] < a2[0]) {return -1;}
																if( a1[0] > a2[0]) {return 1;}
																return 0;
														}});
			
			if( exonVector.firstElement()[2] ==1) {
				
				if( utr5.containsKey(transcript)) {
					
					while(utr5.get(transcript)[1]>=exonVector.firstElement()[1]  ) {
						exonVector.remove(0);
					}
					exonVector.firstElement()[0] = utr5.get(transcript)[1];
					
				}
				if( utr3.containsKey(transcript)) {
					
					while(utr3.get(transcript)[0] <= exonVector.lastElement()[0]) {
						exonVector.remove(exonVector.size()-1);
					}
					
					exonVector.lastElement()[1] = utr3.get(transcript)[0];
				}
				
				int offset = exonVector.firstElement()[0];
				
				for( int i = 0; i< exonVector.size(); i++) {
					exonVector.get(i)[0] = exonVector.get(i)[0] -offset +1;
					exonVector.get(i)[1] = exonVector.get(i)[1] -offset +1;
				}
				
			}else {
				
				
				if( utr5.containsKey(transcript)) {
					while(utr5.get(transcript)[0]<=exonVector.lastElement()[0] ) {
						exonVector.remove(exonVector.size()-1);
					}
					exonVector.lastElement()[1]=utr5.get(transcript)[0];
				}
				
				if( utr3.containsKey(transcript)) {
					while(utr3.get(transcript)[1] >= exonVector.firstElement()[1]) {
						exonVector.remove(0);
					}
					
					exonVector.firstElement()[0] = utr3.get(transcript)[1];
				}
				int offset = exonVector.lastElement()[1];
				Vector<int[]> v = new Vector<int[]>();
				for( int i = exonVector.size()-1; i>=0;i--) {
					int[] a = { offset +1 - exonVector.get(i)[1], offset+1 - exonVector.get(i)[0] };
					v.add(a);
				}
				exonVector = v;
			}
			
			
			Vector<String> motifVector =  motifs.get(transcript);
			
			
			
			
			
			
			String s = transcript + "," + (exonVector.lastElement()[1] ) ;
			
			
			
			
			for(Iterator<String> iterator = motifVector.iterator(); iterator.hasNext();) {
				String motifString = iterator.next();
				String[] split = motifString.split("\\|");
				int[] motif = { (Integer.parseInt(split[1])-1) * 3 +1, (Integer.parseInt(split[2])-1) * 3+3};
				
				int motiflength = motif[1] - motif[0];
					
				int pos = motif[0];
				int toAdd = 0;
				int sumExon = exonVector.get(0)[1] + 1 - exonVector.get(0)[0];
				int index = 0;
				
				
				
				while(pos > sumExon && index < exonVector.size()-1) {
					//System.out.println(transcript + "\t" + split[4] + "\t" + pos  + "\t" + sumExon + "\t" + toAdd + "\t" + exonVector.get(index)[0] + "\t" + exonVector.get(index)[1] + "\t" + exonVector.size());
					toAdd=toAdd + exonVector.get(index+1)[0]- exonVector.get(index)[1] ; 
					sumExon = sumExon + exonVector.get(index)[1] + 1 - exonVector.get(index)[0];
					index++;
				}
				motif[0] = pos + toAdd;
				motif[1] = motif[0] + motiflength;
 				s = s + ",RE|" + motif[0] + "|" + motif[1] + "|" + split[3] + "|" + split[4];
				
			}
			for( int i = 0; i< exonVector.size()-1; i++) {
				int[] intron = {exonVector.get(i)[1], exonVector.get(i+1)[0]  };
				s = s + "," + "RE|" + intron[0] + "|" + intron[1] + "|#ff0000|INTRON";
			}
			out.write(s);
			out.newLine();
		}
		
		out.close();
		
	}
	
	
	public static void writeIntronExonForITOL2(HashSet<String> transcriptsList, HashMap<String,Vector<String>> motifs, File gffHC, File gffLC, File outputFile )throws IOException{
		
		HashMap<String, Transcript> transcripts = new HashMap<String, Transcript>();
		
		File[] gffs = {gffHC, gffLC};
		for( int i = 0; i< gffs.length; i++) {
			BufferedReader in = new BufferedReader(new FileReader(gffs[i]));

			for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
				if(inputline.startsWith("#")) {
					continue;
				}
				String[] split = inputline.split("\t");
				//chr1A	IWGSC_v1.1_201706	exon	40098	40731	.	-	.	ID=TraesCS1A02G000100.1.exon1;Parent=TraesCS1A02G000100.1
				
				if( split[2].equalsIgnoreCase("mrna")) {
					String id = split[8].split("ID=")[1].split(";")[0];
					if(transcriptsList.contains(id)) {
						transcripts.put( id, new Transcript(id) );
						
						if( split[6].equalsIgnoreCase("-")) {
							transcripts.get(id).setForwardStrand(false);
						}else {
							transcripts.get(id).setForwardStrand(true);
						}
					
					}
				}
				if( split[2].equalsIgnoreCase("exon")) {
					String id = split[8].split("Parent=")[1].split(";")[0];
					if(transcriptsList.contains(id)) {
						int[] exon = {Integer.parseInt(split[3]), Integer.parseInt(split[4])};
						transcripts.get(id).addExon(exon);
					}
					
					
				}
				if( split[2].equalsIgnoreCase("cds")) {
					String id = split[8].split("Parent=")[1].split(";")[0];
					if(transcriptsList.contains(id)) {
						int[] cds = {Integer.parseInt(split[3]), Integer.parseInt(split[4])};
						transcripts.get(id).updateCDS(cds);
					}
					
				}
			}
		}
		
		for( Iterator<String> iterator = motifs.keySet().iterator(); iterator.hasNext();) {
			String key = iterator.next();
			if( transcripts.containsKey(key)) {
				transcripts.get(key).addMotifs(motifs.get(key));
			}
		}
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		
		out.write("DATASET_DOMAINS");out.newLine();
		out.write("SEPARATOR COMMA");out.newLine();
		out.write("DATASET_LABEL,Introns + Domains");out.newLine();
		out.write("COLOR,#00ff00");out.newLine();
		out.write("DATASET_SCALE,1000-1000-#0000ff,2000-2000-#ff0000,3000-1000-#00ff00");out.newLine();
		out.write("DATA");out.newLine();
		
		
		
		
		for(Iterator<String> iterator = transcripts.keySet().iterator(); iterator.hasNext();) {
			String key = iterator.next();
			out.write(transcripts.get(key).getItolIntronExonString());
			out.newLine();
		}
		
		out.close();
		
	
	}
	
	
	
	
}	

	

