package readSequenceToRefToMethylation;

import java.io.File;
import java.util.*;

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
				+ "5. Read is first in pair (true/false)\n"
				+ "6. Read mapes to forward strand (true/false)"
				+ "\n\n"
				+ "The strings for read, reference and methylation have the same length"
				+ "\n\n"
				+ "Reads are represented in the same way as they appear on the original FASTQ file\n"
				+ "and the reference is reverse complemented accordingly."
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
		FastaReference fa= new FastaReference( new File(ref) );
		/* ----------------------------------------------------- */
		
		int n= 0;
		for(SAMRecord rec : sam){
			
			// Note use of 'n' to mark soft clipped bases. Hopefully 'n' doesn't appear
			// in read sequence.
			AlignedRead aln= new AlignedRead(rec, fa, (byte)'n'); 
		
			List<Integer> softClippedPositions= aln.getSoftClippedPositions(aln.getReadbases()); 
			byte[] readbases= aln.softClipArray(softClippedPositions, aln.getReadbases());
			byte[] refbases= aln.softClipArray(softClippedPositions, aln.getRefbases());
			byte[] methylbases= aln.softClipArray(softClippedPositions, aln.getMethylbases());
			
			if (readbases.length != refbases.length || refbases.length != methylbases.length){
				System.err.println("Unequal length of read, reference or methylation arrays");
				System.exit(1);
			}
			String readstring= new String(readbases);
			String refstring= new String(refbases);
			String methylstring= new String(methylbases);
			if (readbases.length == 0){
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
			sb.append(aln.isFirst()); 			
			sb.append("\t");
			sb.append(aln.isForward());
			System.out.println( sb.toString() );
						
			n++;
			if (n >= 1000000){
				break;
			}
		}
		sam.close();
	}

}
