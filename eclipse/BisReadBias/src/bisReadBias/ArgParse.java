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
				+ "\nNB: Fasta file is read in memory, adjust memory requirements via -Xmx accordingly"
				+ "\n"
				+ "\nOUTPUT\n"
				+ "To stdout a table with colums methylation counts aggregated along reads"
				+ "\n"
				+ "read position cnt_met cnt_unmet cnt_mismatch cnt_noncytosine\n"
				+ "1    1        5923    29        32           9114\n"
				+ "1    2        4100    12        10           10976\n"
				+ "\nTODO:\n"
				+ "Allow reading bam file from stdin.");	

		parser.addArgument("-f", "--fasta")
			.type(String.class)
			.required(true)
			.help("Input fasta file.");
		parser.addArgument("-i", "--input")
			.type(String.class)
			.required(true)
			.help("Input sam or bam file.");
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
