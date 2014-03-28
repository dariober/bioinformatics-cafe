package fixBSseqSAM;

import java.io.File;
import java.io.IOException;

import net.sf.samtools.*;
import net.sf.samtools.SAMFileReader.ValidationStringency;

public class Main {
	
	
	public static void main(String[] args){

		String bam= args[0];
		String outSam= args[1];
		
		// INPUT
		SAMFileReader sam= bam.equals("-") ? 
				new SAMFileReader ( System.in ) :
				new SAMFileReader ( new  File( bam ));
		sam.setValidationStringency(ValidationStringency.LENIENT);
		
		/* ------------------------------------------------------------------ */

		// Prepare new sam file header
		// ---------------------------
		SAMFileHeader newHeader= null;
		SAMSequenceDictionary newSeqDict= new SAMSequenceDictionary();
		try{
			newHeader= Tools.copyAndGetFileHeader(sam.getFileHeader());
		}
		catch(IOException e) {
			e.printStackTrace();
		}
		newSeqDict= Tools.prepareSequenceDictionary(newHeader.getSequenceDictionary());
		newHeader.setSequenceDictionary(newSeqDict);

		// OUTPUT
		SAMFileWriterFactory sf= new SAMFileWriterFactory();
		SAMFileWriter sw= sf.makeSAMWriter(newHeader, false, System.out );
		if ( !outSam.equals("-") ){
			sw= sf.makeSAMOrBAMWriter(newHeader, false, new File(outSam) );
		}
		
		int n= 0;
		for (SAMRecord rec : sam){
			n += 1;
			rec.setAttribute("YG", Tools.getGenomeConversionTag(rec));
			// Correct reference name
			String newRefName= rec.getReferenceName().replaceAll("_CT$|_GA$", "");
			rec.setReferenceName(newRefName);
			// Correct mate reference name
			String mateReferenceName= rec.getMateReferenceName().replaceAll("_CT$|_GA$", "");
			rec.setMateReferenceName(mateReferenceName); 
			//System.out.println(rec);
			sw.addAlignment(rec);
		}		
		System.err.println(n);
		sam.close();
		sw.close();
		
	}
	
}
