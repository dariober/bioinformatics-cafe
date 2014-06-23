package bisReadBias;

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
				.newArgumentParser("BisReadBias")
				.defaultHelp(true)
				.version("${prog} " + VERSION)
				.description("DESCRIPTION\n" 
				+ "Methylation profile along aligned reads from BS-Seq experiment\n"
				+ "\nUSAGE\n"
				+ "java -Xmx4g BisReadBias.jar -f <ref.fasta> -i <aln.sam>"
				+ "\n"
				+ "samtools view -u <aln.bam> | java -Xmx4g BisReadBias.jar -f <ref.fasta> -i -"
				+ "\n"
				+ "\nNB: Fasta file is read in memory, adjust memory requirements via -Xmx accordingly"
				+ "\n"
				+ "\nOUTPUT\n"
				+ "To stdout a table with colums methylation counts aggregated along reads"
				+ "\nColumns are:"
				+ "\n read            [1|2],"
				+ "\n position        [1-based pos on read],"
				+ "\n cnt_met         [count methylated],"
				+ "\n cnt_unmet       [count unmenthylated],"
				+ "\n cnt_mismatch    [non C or T],"
				+ "\n cnt_noncytosine [not a C in the reference]");	

		parser.addArgument("-f", "--fasta")
			.type(String.class)
			.required(true)
			.help("Input fasta file.");
		parser.addArgument("-i", "--input")
			.type(String.class)
			.required(true)
			.help("Input sam or bam file. Use - to read *bam* from stdin.");
		parser.addArgument("-s", "--step")
			.type(Integer.class)
			.required(false)
			.setDefault(-1)
			.help("Process every s-th alignment from input sam. "
					+ "Default is to process all the reads.");
		
		parser.addArgument("-S", "--stopAfter")
		.type(Integer.class)
		.required(false)
		.setDefault(-1)
		.help("Stop after having collected ths many alignment records. "
				+ "Default is to read all of them.");
		
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
		
		// Optionally used to check args using syntax like:
		// if (opts.getInt("input") == null){
		//		//
		// }

	}
}
