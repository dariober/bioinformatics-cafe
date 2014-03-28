package fixBSseqSAM;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import net.sf.samtools.*;
import net.sf.samtools.SAMFileReader.ValidationStringency;

public class Tools {

	public static SAMFileHeader copyAndGetFileHeader(SAMFileHeader header) throws IOException{

		File temp = File.createTempFile("FixBSseqSAM-FileHeader", ".tmp.sam");
		temp.deleteOnExit();
		
		SAMFileWriterFactory sf= new SAMFileWriterFactory();
		SAMFileWriter sw= sf.makeSAMWriter(header, false, temp );
		sw.close();

		SAMFileReader sam= new SAMFileReader (temp);
		sam.setValidationStringency(ValidationStringency.LENIENT);
		SAMFileHeader newHeader= sam.getFileHeader();
		sam.close();
		return(newHeader);
	}
	
	/**
	 * Generate a new sequence disctionary using the input one as template. The input dict will have 
	 * suffixes stripped from the reference names. 
	 * @param samSeqDict
	 * @return
	 */
	public static SAMSequenceDictionary prepareSequenceDictionary(SAMSequenceDictionary samSeqDict){
		
		List<SAMSequenceRecord> seqList= samSeqDict.getSequences();
		
		// This list is used to know which names have been already put in the new list of sequence names
		List<String> sequenceNames= new ArrayList<String>(); 
		
		List<SAMSequenceRecord> newSeqList= new ArrayList<SAMSequenceRecord>();
		int refIndex= 0; // Reference index to assign to the new names
		for (SAMSequenceRecord rec : seqList){
			
			String newSeqName= rec.getSequenceName().replaceAll("_CT$|_GA$", "");

			if (!sequenceNames.contains(newSeqName)){
				SAMSequenceRecord newRec= new SAMSequenceRecord(newSeqName, rec.getSequenceLength());
				newRec.setSequenceIndex(refIndex);
				newRec.setAssembly(rec.getAssembly());
				newSeqList.add(newRec);
				sequenceNames.add(newSeqName);
				refIndex += 1;
			}
		}		
		SAMSequenceDictionary newSeqDict= new SAMSequenceDictionary(newSeqList);
		return(newSeqDict);
	} 
	
	public static String getGenomeConversionTag(SAMRecord rec){
		String yg= "NA";
		String refName= rec.getReferenceName();
		if (refName.endsWith("_CT")){
			yg= "CT";
		}
		else if(refName.endsWith("_GA")) {
			yg= "GA";
		}
		else if(refName.equals("*")){
			// Leave yg tag as default from above
		}
		else {
			System.err.println("No genome conversion tag found for record " + rec + "; Ref. name " + rec.getReferenceName());
			System.exit(1);
		}
		return(yg);
	}
	
}
