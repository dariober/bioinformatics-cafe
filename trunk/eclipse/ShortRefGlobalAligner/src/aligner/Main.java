package aligner;

import java.io.*;
import java.util.*;
//import java.util.concurrent.ExecutorService;
//import java.util.concurrent.Executors;

import net.sourceforge.argparse4j.inf.Namespace;

import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.DNASequenceCreator;
import org.biojava3.core.sequence.io.FastaReader;
import org.biojava3.core.sequence.io.GenericFastaHeaderParser;

public class Main {

	public static void main(String[] args) throws Exception {

		/* Start parsing arguments */
		Namespace opts= ArgParse.argParse(args);
		ArgParse.validateArgs(opts);
		String fastq= opts.getString("fastq");
		String reference= opts.getString("reference");
		String matrix= opts.getString("matrix");
		String refName= opts.getString("refName");
		int stopAfter= opts.getInt("stopAfter");
		int step= opts.getInt("step");
		
		Aligner aln= new Aligner();
		String matrixFile= null;
		if (matrix.equals("nuc_N0.txt") || matrix.equals("nuc_N1.txt")){
			matrixFile= aln.getTmpFileFromResourceFile(matrix);
		} else {
			// TODO: Code to read custom matrix
		}
		/* End parsing arguments */

		/* Parallelization variables*/
//		int NUM_OF_THREADS= Runtime.getRuntime().availableProcessors();
//		int BATCH_SIZE= 100; // Accumulate this many items (reads) before passing them to parallelize
//      List<String[]> fqreadList= new ArrayList<String[]>(); // Accumulate here reads to pass to parallelize()
//      List<String> outList= new ArrayList<String>();

		
		/* Get input */
		FileInputStream inStream = new FileInputStream( reference );
	
		GenericFastaHeaderParser<DNASequence, NucleotideCompound> genericFastaHeaderParser= new GenericFastaHeaderParser<DNASequence, NucleotideCompound>();
		DNASequenceCreator dnaSequenceCreator= new DNASequenceCreator(AmbiguityDNACompoundSet.getDNACompoundSet());

		FastaReader<DNASequence, NucleotideCompound> fastaReader= new FastaReader<DNASequence, NucleotideCompound>(
			inStream, genericFastaHeaderParser, dnaSequenceCreator); 
		LinkedHashMap<String, DNASequence> b = fastaReader.process();
        String refseq= b.get(refName).getSequenceAsString();

        BufferedReader br= FastqReader.openFastq(fastq);
        String[] fqread= FastqReader.getNextRead(br);
        int nreads= 0;
        int nprocessed= 0;

        while (fqread != null){
        	nreads++;
        	if(nreads % step == 0){
        		String outline= processRead(fqread, refseq, matrixFile);
	        	nprocessed++;
	        	System.out.println(outline);
        	}
        	if(nprocessed >= stopAfter && stopAfter > 0){
        		break;
        	} 
        	fqread= FastqReader.getNextRead(br);
        	if (nreads % 10000 == 0){
        		System.err.println(nreads + " reads");
        	}
		}
        br.close();
        System.err.println(nprocessed + " processed.");
        System.exit(0);
	}
	
	private static String processRead(String[] fqread, String refseq, String matrixFile) throws FileNotFoundException{
		Aligner aln= new Aligner();
		String sequence= fqread[1];
    	aln.align(sequence, refseq, matrixFile);
    	
        FKMod fkmod= new FKMod();
    	String firstBarcode= fkmod.getQueryBasesAtTargetIndexes(aln, FKMod.INDEX_FIRST_BARCODE);
    	String secondBarcode= fkmod.getQueryBasesAtTargetIndexes(aln, FKMod.INDEX_SECOND_BARCODE);
    	String mods= fkmod.getQueryBasesAtTargetIndexes(aln, FKMod.INDEX_MODIFICATIONS);

    	StringBuilder sb= new StringBuilder();
    	sb.append(aln.getPair().getQuery());   sb.append("\t");
    	sb.append(aln.getPair().getTarget());  sb.append("\t");
    	sb.append(aln.getScore());             sb.append("\t");
    	sb.append(aln.getStrand()); 		   sb.append("\t");
    	sb.append(aln.getPair().getLength());  sb.append("\t");
    	sb.append(aln.getPair().getNumIdenticals());  sb.append("\t");
    	sb.append(aln.getPair().getNumSimilars());  sb.append("\t");
    	sb.append(firstBarcode);               sb.append("\t");
    	sb.append(secondBarcode);              sb.append("\t");
    	for(char c : mods.toCharArray()){
    		sb.append(c);             sb.append("\t");
    	}
    	return(sb.toString().trim());
	}
	// See http://stackoverflow.com/questions/4010185/parallel-for-for-java
	// But also
	// See http://stackoverflow.com/questions/5686200/parallelizing-a-for-loop
	/*	
 	private static List<String> parallelize(List<String[]> fqreadList, final String refseq, final String matrixFile, int NUM_OF_THREADS){
		final List<String> outlist= new ArrayList<String>();
		ExecutorService exec = Executors.newFixedThreadPool(NUM_OF_THREADS);
		try {
		    for (final String[] fqread : fqreadList) {
		        exec.submit(new Runnable() {
		            @Override
		            public void run() {
		                // do stuff with o.
						try {
							System.out.println(fqread);
							String outline= processRead(fqread, refseq, matrixFile);							
							outlist.add(outline);
						} catch (FileNotFoundException e) {
							e.printStackTrace();
						}
		            }
		        });
		    }
		} finally {
		    exec.shutdown();
		}
		return outlist;
	}
    */
}
