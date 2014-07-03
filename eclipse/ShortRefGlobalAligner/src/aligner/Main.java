package aligner;

import java.io.*;
import java.util.*;

import net.sourceforge.argparse4j.inf.Namespace;

import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.DNASequenceCreator;
import org.biojava3.core.sequence.io.FastaReader;
import org.biojava3.core.sequence.io.GenericFastaHeaderParser;

// import aligner.ArgParse;

public class Main {

	public static void main(String[] args) throws Exception {

		/* Start parsing arguments */
		Namespace opts= ArgParse.argParse(args);
		ArgParse.validateArgs(opts);
		String fastq= opts.getString("fastq");
		String reference= opts.getString("reference");
		String matrix= opts.getString("matrix");
		String refName= opts.getString("refName");
		
		Aligner aln= new Aligner();
		String matrixFile= null;
		if (matrix.equals("nuc_N0.txt") || matrix.equals("nuc_N1.txt")){
			matrixFile= aln.getTmpFileFromResourceFile(matrix);
		} else {
			// TODO: Code to read custom matrix
		}
		/* End parsing arguments */

		FileInputStream inStream = new FileInputStream( reference );
	
		GenericFastaHeaderParser<DNASequence, NucleotideCompound> genericFastaHeaderParser= new GenericFastaHeaderParser<DNASequence, NucleotideCompound>();
		DNASequenceCreator dnaSequenceCreator= new DNASequenceCreator(AmbiguityDNACompoundSet.getDNACompoundSet());

		FastaReader<DNASequence, NucleotideCompound> fastaReader= new FastaReader<DNASequence, NucleotideCompound>(
			inStream, genericFastaHeaderParser, dnaSequenceCreator); 
		LinkedHashMap<String, DNASequence> b = fastaReader.process();
        String refseq= b.get(refName).getSequenceAsString();

        FKMod fkmod= new FKMod();
        BufferedReader br= FastqReader.openFastq(fastq);
        String[] fqread= FastqReader.getNextRead(br);
        while (fqread != null){
        	String sequence= fqread[1];
        	aln.align(sequence, refseq, matrixFile);
        	
        	String firstBarcode= fkmod.getQueryBasesAtTargetIndexes(aln, FKMod.INDEX_FIRST_BARCODE);
        	String secondBarcode= fkmod.getQueryBasesAtTargetIndexes(aln, FKMod.INDEX_SECOND_BARCODE);
        	String mods= fkmod.getQueryBasesAtTargetIndexes(aln, FKMod.INDEX_MODIFICATIONS);

        	StringBuilder sb= new StringBuilder();
        	sb.append(aln.getPair().getQuery());   sb.append("\t");
        	sb.append(aln.getPair().getTarget());  sb.append("\t");
        	sb.append(aln.getScore());             sb.append("\t");
        	sb.append(aln.getStrand()); 		   sb.append("\t");
        	sb.append(aln.getPair().getLength());  sb.append("\t");
        	sb.append(firstBarcode);               sb.append("\t");
        	sb.append(secondBarcode);              sb.append("\t");
        	for(char c : mods.toCharArray()){
        		sb.append(c);             sb.append("\t");
        	}
        	System.out.println(sb.toString().trim());
        	fqread= FastqReader.getNextRead(br);
		}
        br.close();
	}
}
