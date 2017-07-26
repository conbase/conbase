import os
main_directory = os.path.dirname(os.path.realpath(__file__))
STYLE_PATH = main_directory + '/Style.css'
SCRIPT_PATH = main_directory + '/Script.js'


class TSV(object):
    def __init__(self, path, sample_names):
        self.path = path
        self.sample_names = sample_names
        self.writer = self.open()

    def open(self):
        writer = open(self.path,'w')
        writer.write("CHROM\tPOS\t")
        for sample_name in self.sample_names:
            writer.write(sample_name + "\t")
        writer.write("\n")
        return writer

    def write_site(self, site):
        self.writer.write(site.CHROM + "\t" + str(site.real_POS()) + "\t")
        for sample_name in self.sample_names:
            self.writer.write(site.samples[sample_name].info + "\t")
        self.writer.write("\n")

    def close(self):
        self.writer.close()

class HTML(object):
    def __init__(self, path, sample_names, stats_params, analyze_params, misc_params):
        self.path = path
        self.sample_names = sample_names
        self.source_code = self.open()
        self.stats_params = stats_params
        self.analyze_params = analyze_params
        self.misc_params = misc_params

    def open(self):
        source_code = ['<html>\n<head>\n<style>\n']
        source_code.append(open(STYLE_PATH).read())
        source_code.append("</style>\n<script>\n")
        source_code.append(open(SCRIPT_PATH).read())
        source_code.append("</script>\n</head>\n<body>")
        source_code.append('<body><div id="button"> <button type="button" id="btn" onclick="change(\'min\')">Zoom out</button> <img height="30px" src="https://s17.postimg.org/bb0oek5i7/conbase.gif" style="float:right;  margin-right:20px"><br/></div>')
        source_code.append('\n<table id="data"><tr><td class="cell"> position</td><td class="cell">BULK</td>\n')
        for sample_name in self.sample_names:
            source_code.append('<td class="cell">' + sample_name + '</td>')
        source_code.append('</tr>')
        return source_code
    
    def write_site(self, site):
        self.source_code.append('<tr>\n')

        self.source_code.append('<td class="info cell">chr' + site.CHROM + ':' + str(site.real_POS()) + '</td>\n') 
        bulk_a1_ratio = float(site.BULK_INFO[site.ALTS['A1']])/site.BULK_INFO['SUM']
        if bulk_a1_ratio > 0:
            self.source_code.append('<td class="bulk cell"> <p class="clicker" onclick="show(' + "'hidden_" +  site.CHROM + "_" + str(site.real_POS()) + "'" + ')" > DP: ' + str(site.BULK_INFO['SUM']) + ' (' + '{0:.2f}'.format(bulk_a1_ratio) + ')</p></td>')  
        else:
            self.source_code.append('<td class="bulk cell"> <p class="clicker" onclick="show(' + "'hidden_" +  site.CHROM + "_" + str(site.real_POS()) + "'" + ')" > DP: ' + str(site.BULK_INFO['SUM']) + '</p></td>')  

        for sample_name in self.sample_names:
            cell_name_type = "name"
            if site.samples[sample_name].info == "X":
                cell_name_type = "X"
            elif site.samples[sample_name].info == "HOMO-R" or site.samples[sample_name].info == "ADO-R" or site.samples[sample_name].info == "HOMO-C3":
                cell_name_type = "0"
            elif site.samples[sample_name].info == "HET" or site.samples[sample_name].info == "ADO-A1" or site.samples[sample_name].info == "HET-C3":
                cell_name_type = "1"
            elif site.samples[sample_name].info == "None" and sum(site.samples[sample_name].AD.values()) == 0:
                cell_name_type = "-"
            elif site.samples[sample_name].info == "None" and sum(site.samples[sample_name].AD.values()) > 0:
                cell_name_type = "not informative"
            elif site.samples[sample_name].info == "unknown":
                cell_name_type = "?"
            if site.samples[sample_name].info == "None" and sum(site.samples[sample_name].AD.values()) == 0:
                self.source_code.append('<td class="'+ site.samples[sample_name].info  + ' cell" ><p class="clicker" style="opacity:0.0" onclick="show(' + "'hidden_" +  site.CHROM + "_" + str(site.real_POS()) + "'" + ')" >' + cell_name_type + '</p>\n')
            
            else:
                self.source_code.append('<td class="'+ site.samples[sample_name].info  + ' cell" ><p class="clicker " onclick="show(' + "'hidden_" +  site.CHROM + "_" + str(site.real_POS()) + "'" + ')" >' + cell_name_type + '</p>\n')
            self.source_code.append('<table class="hidden ' + 'hidden_' +  site.CHROM + '_' + str(site.real_POS()) + '">\n')
            self.source_code.append('<tr><td colspan=5><strong>' + sample_name + '</strong>: ' + site.samples[sample_name].get_AD(site) + '</td></tr>')

            self.source_code.append('<tr><td colspan=5>-----------</td></tr>\n')
            msp_list = sorted(list(site.samples[sample_name].MSP.items()), key=lambda x: x[0])
            for snp_pos, msp in msp_list:
                if msp.get_ms_total() > 0:
                    if msp.voted != '':
                        if msp.voted == 'unknown':
                            self.source_code.append('<tr><td class="voted_unknown" colspan=5>SNP: ' + str(snp_pos+1) + ' (ms: ' + str(site.snp_ms_win[snp_pos]) + ', voted: ' + msp.voted + ')</td></tr>\n')
                        else:
                            self.source_code.append('<tr><td class="voted" colspan=5>SNP: ' + str(snp_pos+1) + ' (ms: ' + str(site.snp_ms_win[snp_pos]) + ', voted: ' + msp.voted + ')</td></tr>\n')
                    else:
                        self.source_code.append('<tr><td colspan=5>SNP: ' + str(snp_pos+1) + '</td></tr>\n')

                    self.source_code.append('<tr><td colspan=2>GT distribution</td>\n')

                    self.source_code.append('<td>' + str(msp.homo_R) + '</td>\n')
                    self.source_code.append('<td>' + str(msp.het) + '</td>\n')
                    self.source_code.append('<td>' + str(msp.homo_A1) + '</td></tr>\n')
                    self.source_code.append('<tr><td>MSP</td>\n')
                    self.source_code.append('<td> RR: ' + str(msp.RR) + '</td>\n')
                    self.source_code.append('<td> RA: ' + str(msp.RA) + '</td>\n')
                    self.source_code.append('<td> AR: ' + str(msp.AR) + '</td>\n')
                    self.source_code.append('<td> AA: ' + str(msp.AA) + '</td></tr>\n')
                    self.source_code.append('<tr><td colspan=5>-----------</td></tr>\n')
            self.source_code.append('</table>\n')
            self.source_code.append('</td>')
        self.source_code.append('</tr>\n')
       
        
    def write_params(self):
        self.source_code.append('\n<br><br><br><table class="params">\n')
        rows = max(len(self.stats_params), len(self.analyze_params), len(self.misc_params))
        stats = list(self.stats_params.items())
        analyze = list(self.analyze_params.items())
        misc = list(self.misc_params.items())
        self.source_code.append('\n<tr class="params_name_row"><td></td><td>Stats params</td><td></td><td>Analyze params</td><td></td><td>Misc params</td><td></td></tr>\n')

        for i in range(rows):
            self.source_code.append('\n<tr><td></td>\n')
            if stats:
                key, value = stats.pop(0)
                self.source_code.append('<td><p>' + str(key) + ':</p></td><td><p>' + str(value) + '</p></td>')
            else:
                self.source_code.append('<td></td><td></td>')
            if analyze:
                key, value = analyze.pop(0)
                self.source_code.append('<td><p>' + str(key) + ':</p></td><td><p>' + str(value) + '</p></td>')
            else:
                self.source_code.append('<td></td><td></td>')
            if misc:
                key, value = misc.pop(0)
                self.source_code.append('<td><p>' + str(key) + ':</p></td><td><p>' + str(value) + '</p></td>')
            else:
                self.source_code.append('<td></td><td></td>')
            self.source_code.append('\n</tr>\n')    
        self.source_code.append('\n</table>\n')

    def close(self):
        self.source_code.append('\n</table>\n')
        self.write_params()
        self.source_code.append('</body>\n</html>')
        f = open(self.path,'w')
        f.write("".join(self.source_code))
        f.close()
        
class Site2JSON(object):
    def __init__(self, path):
        self.f = open(path,'w')    
    def write(self,site):
        if site.TYPE == '':
            self.f.write(repr(site) + '\n')

    def close(self):
        self.f.close()
