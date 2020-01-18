import os


class GenomeVersion(object):

    def __init__(self, config):
        self.config = config

    @property
    def sequence_filename(self):
        result = os.path.join(
            self.config.genome.work_dir,
            'genome.fa'
        )
        assert os.path.exists(result), result
        return result

    @property
    def index_prefix(self):
        result = os.path.join(
            self.config.genome.work_dir,
            self.config.genome.index
        )
        return result

    @staticmethod
    def from_config(config):
        if config.genome.version == 'hg19':
            return HumanGenome19(config)
        elif config.genome.version == 'hg38':
            return HumanGenome38(config)

        raise NotImplementedError()


class HumanGenome19(GenomeVersion):

    def __init__(self, config):
        super(HumanGenome19, self).__init__(config)

    VERSION = "hg19"

    CHROMS = [
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
        "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
        "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
        "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"
    ]
    CHROMS_ALL = [
        'chr1',
        'chr1_gl000191_random',
        'chr1_gl000192_random',
        'chr2',
        'chr3',
        # 'chr4_ctg9_hap1',
        'chr4',
        'chr4_gl000193_random',
        'chr4_gl000194_random',
        'chr5',
        # 'chr6_apd_hap1',
        # 'chr6_cox_hap2',
        # 'chr6_dbb_hap3',
        'chr6',
        # 'chr6_mann_hap4',
        # 'chr6_mcf_hap5',
        # 'chr6_qbl_hap6',
        # 'chr6_ssto_hap7',
        'chr7',
        'chr7_gl000195_random',
        'chr8',
        'chr8_gl000196_random',
        'chr8_gl000197_random',
        'chr9',
        'chr9_gl000198_random',
        'chr9_gl000199_random',
        'chr9_gl000200_random',
        'chr9_gl000201_random',
        'chr10',
        'chr11',
        'chr11_gl000202_random',
        'chr12',
        'chr13',
        'chr14',
        'chr15',
        'chr16',
        # 'chr17_ctg5_hap1',
        'chr17',
        'chr17_gl000203_random',
        'chr17_gl000204_random',
        'chr17_gl000205_random',
        'chr17_gl000206_random',
        'chr18',
        'chr18_gl000207_random',
        'chr19',
        'chr19_gl000208_random',
        'chr19_gl000209_random',
        'chr20',
        'chr21',
        'chr21_gl000210_random',
        'chr22',
        'chrM',
        'chrUn_gl000211',
        'chrUn_gl000212',
        'chrUn_gl000213',
        'chrUn_gl000214',
        'chrUn_gl000215',
        'chrUn_gl000216',
        'chrUn_gl000217',
        'chrUn_gl000218',
        'chrUn_gl000219',
        'chrUn_gl000220',
        'chrUn_gl000221',
        'chrUn_gl000222',
        'chrUn_gl000223',
        'chrUn_gl000224',
        'chrUn_gl000225',
        'chrUn_gl000226',
        'chrUn_gl000227',
        'chrUn_gl000228',
        'chrUn_gl000229',
        'chrUn_gl000230',
        'chrUn_gl000231',
        'chrUn_gl000232',
        'chrUn_gl000233',
        'chrUn_gl000234',
        'chrUn_gl000235',
        'chrUn_gl000236',
        'chrUn_gl000237',
        'chrUn_gl000238',
        'chrUn_gl000239',
        'chrUn_gl000240',
        'chrUn_gl000241',
        'chrUn_gl000242',
        'chrUn_gl000243',
        'chrUn_gl000244',
        'chrUn_gl000245',
        'chrUn_gl000246',
        'chrUn_gl000247',
        'chrUn_gl000248',
        'chrUn_gl000249',
        'chrX',
        'chrY',
    ]

    CHRY_PARs = [
        (10000, 2649520),
        (59034049, 59363566),
    ]


class HumanGenome38(GenomeVersion):

    def __init__(self, config):
        super(HumanGenome19, self).__init__(config)

    VERSION = "hg38"

    CHROMS = [
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
        "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
        "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
        "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"
    ]

    CHRY_PARs = [
        (10000, 2781479),
        (56887902, 57217415),
    ]

    CHROMS_ALL = [
        "chr10",
        "chr10_GL383545v1_alt",
        "chr10_GL383546v1_alt",
        "chr10_KI270824v1_alt",
        "chr10_KI270825v1_alt",
        "chr11",
        "chr11_GL383547v1_alt",
        "chr11_JH159136v1_alt",
        "chr11_JH159137v1_alt",
        "chr11_KI270721v1_random",
        "chr11_KI270826v1_alt",
        "chr11_KI270827v1_alt",
        "chr11_KI270829v1_alt",
        "chr11_KI270830v1_alt",
        "chr11_KI270831v1_alt",
        "chr11_KI270832v1_alt",
        "chr11_KI270902v1_alt",
        "chr11_KI270903v1_alt",
        "chr11_KI270927v1_alt",
        "chr12",
        "chr12_GL383549v1_alt",
        "chr12_GL383550v2_alt",
        "chr12_GL383551v1_alt",
        "chr12_GL383552v1_alt",
        "chr12_GL383553v2_alt",
        "chr12_GL877875v1_alt",
        "chr12_GL877876v1_alt",
        "chr12_KI270833v1_alt",
        "chr12_KI270834v1_alt",
        "chr12_KI270835v1_alt",
        "chr12_KI270836v1_alt",
        "chr12_KI270837v1_alt",
        "chr12_KI270904v1_alt",
        "chr13",
        "chr13_KI270838v1_alt",
        "chr13_KI270839v1_alt",
        "chr13_KI270840v1_alt",
        "chr13_KI270841v1_alt",
        "chr13_KI270842v1_alt",
        "chr13_KI270843v1_alt",
        "chr14",
        "chr14_GL000009v2_random",
        "chr14_GL000194v1_random",
        "chr14_GL000225v1_random",
        "chr14_KI270722v1_random",
        "chr14_KI270723v1_random",
        "chr14_KI270724v1_random",
        "chr14_KI270725v1_random",
        "chr14_KI270726v1_random",
        "chr14_KI270844v1_alt",
        "chr14_KI270845v1_alt",
        "chr14_KI270846v1_alt",
        "chr14_KI270847v1_alt",
        "chr15",
        "chr15_GL383554v1_alt",
        "chr15_GL383555v2_alt",
        "chr15_KI270727v1_random",
        "chr15_KI270848v1_alt",
        "chr15_KI270849v1_alt",
        "chr15_KI270850v1_alt",
        "chr15_KI270851v1_alt",
        "chr15_KI270852v1_alt",
        "chr15_KI270905v1_alt",
        "chr15_KI270906v1_alt",
        "chr16",
        "chr16_GL383556v1_alt",
        "chr16_GL383557v1_alt",
        "chr16_KI270728v1_random",
        "chr16_KI270853v1_alt",
        "chr16_KI270854v1_alt",
        "chr16_KI270855v1_alt",
        "chr16_KI270856v1_alt",
        "chr17",
        "chr17_GL000205v2_random",
        "chr17_GL000258v2_alt",
        "chr17_GL383563v3_alt",
        "chr17_GL383564v2_alt",
        "chr17_GL383565v1_alt",
        "chr17_GL383566v1_alt",
        "chr17_JH159146v1_alt",
        "chr17_JH159147v1_alt",
        "chr17_JH159148v1_alt",
        "chr17_KI270729v1_random",
        "chr17_KI270730v1_random",
        "chr17_KI270857v1_alt",
        "chr17_KI270858v1_alt",
        "chr17_KI270859v1_alt",
        "chr17_KI270860v1_alt",
        "chr17_KI270861v1_alt",
        "chr17_KI270862v1_alt",
        "chr17_KI270907v1_alt",
        "chr17_KI270908v1_alt",
        "chr17_KI270909v1_alt",
        "chr17_KI270910v1_alt",
        "chr18",
        "chr18_GL383567v1_alt",
        "chr18_GL383568v1_alt",
        "chr18_GL383569v1_alt",
        "chr18_GL383570v1_alt",
        "chr18_GL383571v1_alt",
        "chr18_GL383572v1_alt",
        "chr18_KI270863v1_alt",
        "chr18_KI270864v1_alt",
        "chr18_KI270911v1_alt",
        "chr18_KI270912v1_alt",
        "chr19",
        "chr19_GL000209v2_alt",
        "chr19_GL383573v1_alt",
        "chr19_GL383574v1_alt",
        "chr19_GL383575v2_alt",
        "chr19_GL383576v1_alt",
        "chr19_GL949746v1_alt",
        "chr19_GL949747v2_alt",
        "chr19_GL949748v2_alt",
        "chr19_GL949749v2_alt",
        "chr19_GL949750v2_alt",
        "chr19_GL949751v2_alt",
        "chr19_GL949752v1_alt",
        "chr19_GL949753v2_alt",
        "chr19_KI270865v1_alt",
        "chr19_KI270866v1_alt",
        "chr19_KI270867v1_alt",
        "chr19_KI270868v1_alt",
        "chr19_KI270882v1_alt",
        "chr19_KI270883v1_alt",
        "chr19_KI270884v1_alt",
        "chr19_KI270885v1_alt",
        "chr19_KI270886v1_alt",
        "chr19_KI270887v1_alt",
        "chr19_KI270888v1_alt",
        "chr19_KI270889v1_alt",
        "chr19_KI270890v1_alt",
        "chr19_KI270891v1_alt",
        "chr19_KI270914v1_alt",
        "chr19_KI270915v1_alt",
        "chr19_KI270916v1_alt",
        "chr19_KI270917v1_alt",
        "chr19_KI270918v1_alt",
        "chr19_KI270919v1_alt",
        "chr19_KI270920v1_alt",
        "chr19_KI270921v1_alt",
        "chr19_KI270922v1_alt",
        "chr19_KI270923v1_alt",
        "chr19_KI270929v1_alt",
        "chr19_KI270930v1_alt",
        "chr19_KI270931v1_alt",
        "chr19_KI270932v1_alt",
        "chr19_KI270933v1_alt",
        "chr19_KI270938v1_alt",
        "chr1",
        "chr1_GL383518v1_alt",
        "chr1_GL383519v1_alt",
        "chr1_GL383520v2_alt",
        "chr1_KI270706v1_random",
        "chr1_KI270707v1_random",
        "chr1_KI270708v1_random",
        "chr1_KI270709v1_random",
        "chr1_KI270710v1_random",
        "chr1_KI270711v1_random",
        "chr1_KI270712v1_random",
        "chr1_KI270713v1_random",
        "chr1_KI270714v1_random",
        "chr1_KI270759v1_alt",
        "chr1_KI270760v1_alt",
        "chr1_KI270761v1_alt",
        "chr1_KI270762v1_alt",
        "chr1_KI270763v1_alt",
        "chr1_KI270764v1_alt",
        "chr1_KI270765v1_alt",
        "chr1_KI270766v1_alt",
        "chr1_KI270892v1_alt",
        "chr20",
        "chr20_GL383577v2_alt",
        "chr20_KI270869v1_alt",
        "chr20_KI270870v1_alt",
        "chr20_KI270871v1_alt",
        "chr21",
        "chr21_GL383578v2_alt",
        "chr21_GL383579v2_alt",
        "chr21_GL383580v2_alt",
        "chr21_GL383581v2_alt",
        "chr21_KI270872v1_alt",
        "chr21_KI270873v1_alt",
        "chr21_KI270874v1_alt",
        "chr22",
        "chr22_GL383582v2_alt",
        "chr22_GL383583v2_alt",
        "chr22_KB663609v1_alt",
        "chr22_KI270731v1_random",
        "chr22_KI270732v1_random",
        "chr22_KI270733v1_random",
        "chr22_KI270734v1_random",
        "chr22_KI270735v1_random",
        "chr22_KI270736v1_random",
        "chr22_KI270737v1_random",
        "chr22_KI270738v1_random",
        "chr22_KI270739v1_random",
        "chr22_KI270875v1_alt",
        "chr22_KI270876v1_alt",
        "chr22_KI270877v1_alt",
        "chr22_KI270878v1_alt",
        "chr22_KI270879v1_alt",
        "chr22_KI270928v1_alt",
        "chr2",
        "chr2_GL383521v1_alt",
        "chr2_GL383522v1_alt",
        "chr2_GL582966v2_alt",
        "chr2_KI270715v1_random",
        "chr2_KI270716v1_random",
        "chr2_KI270767v1_alt",
        "chr2_KI270768v1_alt",
        "chr2_KI270769v1_alt",
        "chr2_KI270770v1_alt",
        "chr2_KI270771v1_alt",
        "chr2_KI270772v1_alt",
        "chr2_KI270773v1_alt",
        "chr2_KI270774v1_alt",
        "chr2_KI270775v1_alt",
        "chr2_KI270776v1_alt",
        "chr2_KI270893v1_alt",
        "chr2_KI270894v1_alt",
        "chr3",
        "chr3_GL000221v1_random",
        "chr3_GL383526v1_alt",
        "chr3_JH636055v2_alt",
        "chr3_KI270777v1_alt",
        "chr3_KI270778v1_alt",
        "chr3_KI270779v1_alt",
        "chr3_KI270780v1_alt",
        "chr3_KI270781v1_alt",
        "chr3_KI270782v1_alt",
        "chr3_KI270783v1_alt",
        "chr3_KI270784v1_alt",
        "chr3_KI270895v1_alt",
        "chr3_KI270924v1_alt",
        "chr3_KI270934v1_alt",
        "chr3_KI270935v1_alt",
        "chr3_KI270936v1_alt",
        "chr3_KI270937v1_alt",
        "chr4",
        "chr4_GL000008v2_random",
        "chr4_GL000257v2_alt",
        "chr4_GL383527v1_alt",
        "chr4_GL383528v1_alt",
        "chr4_KI270785v1_alt",
        "chr4_KI270786v1_alt",
        "chr4_KI270787v1_alt",
        "chr4_KI270788v1_alt",
        "chr4_KI270789v1_alt",
        "chr4_KI270790v1_alt",
        "chr4_KI270896v1_alt",
        "chr4_KI270925v1_alt",
        "chr5",
        "chr5_GL000208v1_random",
        "chr5_GL339449v2_alt",
        "chr5_GL383530v1_alt",
        "chr5_GL383531v1_alt",
        "chr5_GL383532v1_alt",
        "chr5_GL949742v1_alt",
        "chr5_KI270791v1_alt",
        "chr5_KI270792v1_alt",
        "chr5_KI270793v1_alt",
        "chr5_KI270794v1_alt",
        "chr5_KI270795v1_alt",
        "chr5_KI270796v1_alt",
        "chr5_KI270897v1_alt",
        "chr5_KI270898v1_alt",
        "chr6",
        "chr6_GL000250v2_alt",
        "chr6_GL000251v2_alt",
        "chr6_GL000252v2_alt",
        "chr6_GL000253v2_alt",
        "chr6_GL000254v2_alt",
        "chr6_GL000255v2_alt",
        "chr6_GL000256v2_alt",
        "chr6_GL383533v1_alt",
        "chr6_KB021644v2_alt",
        "chr6_KI270758v1_alt",
        "chr6_KI270797v1_alt",
        "chr6_KI270798v1_alt",
        "chr6_KI270799v1_alt",
        "chr6_KI270800v1_alt",
        "chr6_KI270801v1_alt",
        "chr6_KI270802v1_alt",
        "chr7",
        "chr7_GL383534v2_alt",
        "chr7_KI270803v1_alt",
        "chr7_KI270804v1_alt",
        "chr7_KI270805v1_alt",
        "chr7_KI270806v1_alt",
        "chr7_KI270807v1_alt",
        "chr7_KI270808v1_alt",
        "chr7_KI270809v1_alt",
        "chr7_KI270899v1_alt",
        "chr8",
        "chr8_KI270810v1_alt",
        "chr8_KI270811v1_alt",
        "chr8_KI270812v1_alt",
        "chr8_KI270813v1_alt",
        "chr8_KI270814v1_alt",
        "chr8_KI270815v1_alt",
        "chr8_KI270816v1_alt",
        "chr8_KI270817v1_alt",
        "chr8_KI270818v1_alt",
        "chr8_KI270819v1_alt",
        "chr8_KI270820v1_alt",
        "chr8_KI270821v1_alt",
        "chr8_KI270822v1_alt",
        "chr8_KI270900v1_alt",
        "chr8_KI270901v1_alt",
        "chr8_KI270926v1_alt",
        "chr9",
        "chr9_GL383539v1_alt",
        "chr9_GL383540v1_alt",
        "chr9_GL383541v1_alt",
        "chr9_GL383542v1_alt",
        "chr9_KI270717v1_random",
        "chr9_KI270718v1_random",
        "chr9_KI270719v1_random",
        "chr9_KI270720v1_random",
        "chr9_KI270823v1_alt",
        "chrM",
        "chrUn_GL000195v1",
        "chrUn_GL000213v1",
        "chrUn_GL000214v1",
        "chrUn_GL000216v2",
        "chrUn_GL000218v1",
        "chrUn_GL000219v1",
        "chrUn_GL000220v1",
        "chrUn_GL000224v1",
        "chrUn_GL000226v1",
        "chrUn_KI270302v1",
        "chrUn_KI270303v1",
        "chrUn_KI270304v1",
        "chrUn_KI270305v1",
        "chrUn_KI270310v1",
        "chrUn_KI270311v1",
        "chrUn_KI270312v1",
        "chrUn_KI270315v1",
        "chrUn_KI270316v1",
        "chrUn_KI270317v1",
        "chrUn_KI270320v1",
        "chrUn_KI270322v1",
        "chrUn_KI270329v1",
        "chrUn_KI270330v1",
        "chrUn_KI270333v1",
        "chrUn_KI270334v1",
        "chrUn_KI270335v1",
        "chrUn_KI270336v1",
        "chrUn_KI270337v1",
        "chrUn_KI270338v1",
        "chrUn_KI270340v1",
        "chrUn_KI270362v1",
        "chrUn_KI270363v1",
        "chrUn_KI270364v1",
        "chrUn_KI270366v1",
        "chrUn_KI270371v1",
        "chrUn_KI270372v1",
        "chrUn_KI270373v1",
        "chrUn_KI270374v1",
        "chrUn_KI270375v1",
        "chrUn_KI270376v1",
        "chrUn_KI270378v1",
        "chrUn_KI270379v1",
        "chrUn_KI270381v1",
        "chrUn_KI270382v1",
        "chrUn_KI270383v1",
        "chrUn_KI270384v1",
        "chrUn_KI270385v1",
        "chrUn_KI270386v1",
        "chrUn_KI270387v1",
        "chrUn_KI270388v1",
        "chrUn_KI270389v1",
        "chrUn_KI270390v1",
        "chrUn_KI270391v1",
        "chrUn_KI270392v1",
        "chrUn_KI270393v1",
        "chrUn_KI270394v1",
        "chrUn_KI270395v1",
        "chrUn_KI270396v1",
        "chrUn_KI270411v1",
        "chrUn_KI270412v1",
        "chrUn_KI270414v1",
        "chrUn_KI270417v1",
        "chrUn_KI270418v1",
        "chrUn_KI270419v1",
        "chrUn_KI270420v1",
        "chrUn_KI270422v1",
        "chrUn_KI270423v1",
        "chrUn_KI270424v1",
        "chrUn_KI270425v1",
        "chrUn_KI270429v1",
        "chrUn_KI270435v1",
        "chrUn_KI270438v1",
        "chrUn_KI270442v1",
        "chrUn_KI270448v1",
        "chrUn_KI270465v1",
        "chrUn_KI270466v1",
        "chrUn_KI270467v1",
        "chrUn_KI270468v1",
        "chrUn_KI270507v1",
        "chrUn_KI270508v1",
        "chrUn_KI270509v1",
        "chrUn_KI270510v1",
        "chrUn_KI270511v1",
        "chrUn_KI270512v1",
        "chrUn_KI270515v1",
        "chrUn_KI270516v1",
        "chrUn_KI270517v1",
        "chrUn_KI270518v1",
        "chrUn_KI270519v1",
        "chrUn_KI270521v1",
        "chrUn_KI270522v1",
        "chrUn_KI270528v1",
        "chrUn_KI270529v1",
        "chrUn_KI270530v1",
        "chrUn_KI270538v1",
        "chrUn_KI270539v1",
        "chrUn_KI270544v1",
        "chrUn_KI270548v1",
        "chrUn_KI270579v1",
        "chrUn_KI270580v1",
        "chrUn_KI270581v1",
        "chrUn_KI270582v1",
        "chrUn_KI270583v1",
        "chrUn_KI270584v1",
        "chrUn_KI270587v1",
        "chrUn_KI270588v1",
        "chrUn_KI270589v1",
        "chrUn_KI270590v1",
        "chrUn_KI270591v1",
        "chrUn_KI270593v1",
        "chrUn_KI270741v1",
        "chrUn_KI270742v1",
        "chrUn_KI270743v1",
        "chrUn_KI270744v1",
        "chrUn_KI270745v1",
        "chrUn_KI270746v1",
        "chrUn_KI270747v1",
        "chrUn_KI270748v1",
        "chrUn_KI270749v1",
        "chrUn_KI270750v1",
        "chrUn_KI270751v1",
        "chrUn_KI270752v1",
        "chrUn_KI270753v1",
        "chrUn_KI270754v1",
        "chrUn_KI270755v1",
        "chrUn_KI270756v1",
        "chrUn_KI270757v1",
        "chrX",
        "chrX_KI270880v1_alt",
        "chrX_KI270881v1_alt",
        "chrX_KI270913v1_alt",
        "chrY",
        "chrY_KI270740v1_random",        
    ]

