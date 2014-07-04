package aligner;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;


public class ArgParse {
	
	public static String VERSION= "0.1.0";
	
	/* Parse command line args */
	public static Namespace argParse(String[] args){
		ArgumentParser parser= ArgumentParsers
				.newArgumentParser("ShortRefGlobalAligner")
				.defaultHelp(true)
				.version("${prog} " + VERSION)
				.description("DESCRIPTION\n"
						+ "Align short queries to a single reference sequence via global alignment.\n"
						+ "This progrm was writen to align Fumiko's reads against the modifed control.\n"
						+ ""
						+ "Output format:\n"
						+ "1. Aligned query. Original sequence reverse complemented if aligned to reverse strand\n"
						+ "2. Aligned reference\n"
						+ "3. Score\n"
						+ "4. Strand\n"
						+ "5. Aligment length\n"
						+ "6. Num identical\n"
						+ "7. Num similar\n"
						+ "8. Sequence of first barcode\n"
						+ "9. Sequence of second barcode\n"
						+ "10+ Query base at given positions on the reference.\n"
						+ "\n"
						+ "TODO: Args to get custom barcodes and specific positions where modified bases are");
		parser.addArgument("-f", "--fastq")
			.type(String.class)
			.required(true)
			.help("Input fastq file. Can be gzipped.");
		parser.addArgument("-r", "--reference")
			.type(String.class)
			.required(true)
			.help("Reference fasta file. Only one sequence will be used.");
		parser.addArgument("-n", "--refName")
			.type(String.class)
			.setDefault("fk_mod10")
			.help("Name of sequence in fasta file to use for alignment.");
		parser.addArgument("-m", "--matrix")
			.type(String.class)
			.setDefault("nuc_N0.txt")
			.help("File with scoring matrix. Default is read from resource file. Note that in the default matrix N has 0 score.");
		parser.addArgument("-S", "--stopAfter")
			.type(Integer.class)
			.setDefault(-1)
			.help("Stop after having processed this many reads. Default is to process entire fastq file.");
		parser.addArgument("-s", "--step")
			.type(Integer.class)
			.setDefault(1)
			.help("Process every s'th read.");

		
		parser.addArgument("--version", "-v").action(Arguments.version());

		Namespace opts= null;
		try{
			opts= parser.parseArgs(args);
		}
		catch(ArgumentParserException e) {
			parser.handleError(e);
			System.exit(1);
		}		
		return(opts);
	}
	
	public static void validateArgs(Namespace opts){
		//		
	}
}
