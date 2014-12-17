package sequenceMatcher;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;

import org.biojava3.core.sequence.DNASequence;

import net.sourceforge.argparse4j.inf.Namespace;

public class Main {

	public static void main(String[] args) throws IOException{

		/* Start parsing arguments */
		Namespace opts= ArgParse.argParse(args);
		ArgParse.validateArgs(opts);
		
		String subcmd= opts.getString("subcmd");
		
		if(subcmd.equals("match")){
			matcher(opts);
		} else if(subcmd.equals("convert")){
			System.out.println("Not implemented yet");
		} else {
			System.out.println("Unrecognized cmd.");
		}
		
	}

	/**
	 * Routine to convert sam to tab and back. Called with SequenceMatcher convert ...
	 * @param opts Arguments taken from command line and parsed by ArgParse class
	 * @throws IOException 
	 */
	private static void converter(Namespace opts) throws IOException{

		String input= opts.getString("input");
		String a= opts.getString("a");
		String outfmt= opts.getString("outfmt");
		
		if(outfmt.equals("outfmt")){
		// * Convert tab to sam: 
		// ** Read reference file, output as sam header
			System.out.println("@HD\tVN:1.0");
			System.out.println(Sam.fastaFileToSQHeader(a));
			
			// ** Read input tab line by line.
			BufferedReader br= Opener.openBr(input);
			String line;
			while((line= br.readLine()) != null){
				// Skip header if present:
				if(line.startsWith("seq_A\tseq_B\tstrand\t")){
					// Skip header if present:
					continue;
				}
				// ** Convert each line to Match object
				// ** Convert Match object to Sam object
				// ** Output Sam object as sam line.
			}
		}
	}
	
	/**
	 * Routine to perform and output matching. Called with SequenceMatcher match ...
	 * @param opts Arguments taken from command line and parsed by ArgParse class
	 * @throws IOException
	 */
	private static void matcher(Namespace opts) throws IOException{
		
		/* Start parsing arguments */
		String a= opts.getString("a");
		String b= opts.getString("b");
		String method= opts.getString("method");
		int nm= opts.getInt("nm");
		boolean norc= opts.getBoolean("norc");
		boolean noaln= opts.getBoolean("noaln");
		boolean noLD= opts.getBoolean("noLD");
		boolean noJWD= opts.getBoolean("noJWD");
		String outfmt= opts.getString("outfmt");
		
		/* Adjust args */
		
		// How many loops thorough each sequence in B? 
		// 1 if rev comp is not required. 2 otherwise (+ and -)
		int nLoops= (!norc) ? 2 : 1;
		
		// ---------------------------------------------------------------------
		
		ArrayList<String[]> fastaFileA= SequenceReader.readFastaToList(a);
		
		System.err.println("Sequences in A (" + a + "): " + fastaFileA.size());
		
		// ---------------------------------------------------------------------
		// Prepare reading file b
		BufferedReader brB= Opener.openBr(b);
				
		long t0= System.currentTimeMillis();

		/* -------------------------------------------------------------------- */
		// Prepare header:
		if(outfmt.equals("sam")){
			System.out.println("@HD\tVN:1.0");
			System.out.println(Sam.fastaListToSQHeader(fastaFileA));
		} else if(outfmt.equals("tab")){
			System.out.println(new Match().getHeader());
		} else {
			System.err.println("Unsupproted output format");
			System.exit(1);
		}
		/* -------------------------------------------------------------------- */
		int nseq= 0;
		int nmatch= 0;
		String[] seqB;
		String seqBrc= null;
		while((seqB= SequenceReader.getNextSequence(brB)) != null){

			if(!norc){
				seqBrc= new DNASequence(seqB[1]).getReverseComplement().getSequenceAsString();
			}
			
			for(int i=0; i < fastaFileA.size(); i++){
				int nLoopCnt= 0;
								
				while(nLoopCnt < nLoops) {
					String[] seqA= fastaFileA.get(i);
					Match m= new Match(seqA, seqB);
					
					if(nLoopCnt == 0){
						m.setStrand("+");
					} else if (nLoopCnt == 1){
						m.setStrand("-");
						m.setSeqB(seqBrc);
					} else {
						System.exit(1);
					}
					double d= m.getFilterDistance(method, nm);
					
					if(nm < 0 || (d >= 0 && d <= nm)){
						// Edit distance is below threshold. 
						// Fill up remaining metrics:
						if(method != "LD" && !noLD){
							m.computeLD();
						} 
						if(method != "HD"){
							m.computeHD();
						}
						if(!noJWD){
							m.computeJWD();
						}
						if(!noaln){
							m.align();
						}
						nmatch++;
						if(outfmt.equals("sam")){
							System.out.println(Sam.matchToSam(m).toString());
						}else{
							System.out.println(m);
						}
					}
					nLoopCnt++;
				}
			} // end for loop file A
			nseq++;
			if(nseq % 1000 == 0){
				System.err.println("Processed " + nseq + " sequences from B");
			}
		}
		
		long t1= System.currentTimeMillis();
		System.err.println("Completed in " + (t1-t0)/1000.0 + "s");
		System.err.println("N matches: " + nmatch);		
		System.exit(0);
		
	}
	
}
