package readSequenceToRefToMethylation;

import java.io.File;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;

public class Main {

	public static void main(String[] args) {
		
		String docstring= ""
				+ "DESCRIPTION"
				+ "\n"
				+ "For each read in a bam file, returns a tab separated row with:"
				+ "\n"
				+ "1. Read name\n"
				+ "2. String of read bases\n"
				+ "3. String of reference bases for each read base\n"
				+ "4. String representing the methylation status of each base in the read.\n"
				+ "   M= Methylated; u= unmethylated; .= Not a C; *= Mismatch (not C or T)\n"
				+ "5. Read maps to forward strand (true/false)"
				+ "6. Read is first in pair (true/false)\n"
				+ "\n\n"
				+ "The strings for read, reference and methylation have the same length"
				+ "\n\n"
				+ "Soft clipped bases are trimmed and if a read is completely clipped is represented as '!'"
				+ "\n"
				+ "See AlignedRead.java for up do date coding of methylation status"
				+ "\n\n"
				+ "USAGE"
				+ "\n"
				+ "ReadSequenceToRefToMethylation.jar <BAM|SAM> <reference-fasta>"
				+ "\n\n"
				+ "Use - to read alignment from stdin."
				+ "\n\n"
				+ "Version: 0.1.0";
		
		if (args.length != 2 || args[0].equals("-h") || args[0].equals("--help")){
			System.err.println(docstring);
			System.exit(1);
		}
		
		String bam= args[0];
		String ref= args[1];
		/* ----------------------------------------------------- */
		
		// INPUT BAM
		SAMFileReader sam= bam.equals("-") ? 
				new SAMFileReader ( System.in ) :
				new SAMFileReader ( new  File( bam ));
		sam.setValidationStringency(ValidationStringency.SILENT);

		// INPUT FASTA
		System.err.print("Reading reference '" + ref + "'... ");
		FastaReference fa= new FastaReference( new File(ref) );
		System.err.println("Done");
		/* ----------------------------------------------------- */
		
		int n= 0;
		for(SAMRecord rec : sam){
			
			// Note use of 'n' to mark soft clipped bases. Hopefully 'n' doesn't appear
			// in read sequence.
			AlignedRead aln= new AlignedRead(rec, fa, (byte)'n'); 
			aln.clipAlignment();

			if (aln.getReadbases().length != aln.getRefbases().length || aln.getRefbases().length != aln.getMethylbases().length){
				System.err.println("Unequal length of read, reference or methylation arrays");
				System.exit(1);
			}
			String readstring= new String(aln.getReadbases());
			String refstring= new String(aln.getRefbases());
			String methylstring= new String(aln.getMethylbases());
			if (aln.getReadbases().length == 0){
				readstring= "!"; // This is to mark empty reads (completely sof clipped). Use a string that cannot be found in read/ref or methyl (i.e. '*' or 'NA' are NOT suitable)
				refstring= "!";
				methylstring= "!";
			}
			
			StringBuilder sb= new StringBuilder();
			sb.append(rec.getReadName()); 		
			sb.append("\t");
			sb.append(readstring); 	
			sb.append("\t");
			sb.append(refstring); 		
			sb.append("\t");
			sb.append(methylstring); 
			sb.append("\t");
			sb.append(aln.isForward());
			sb.append("\t");
			sb.append(aln.isFirst()); 			
			System.out.println( sb.toString() );
						
			n++;
			if (n % 1000000 == 0){
				System.err.println(n + " reads processed");
			}
		}
		sam.close();
		System.err.println(n + " reads processed");
	}

}
