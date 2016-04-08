csx_server_path = "http://consensx.itk.ppke.hu/public/_include/"

csx_server_path += '/'

def writeHeaderHTML(path, version):
    html = path + "result_sheet.php"
    html = open(html, 'w')

    html.write("""
<!DOCTYPE html>
<html class="" lang="en">
<head>
  <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
  <title>CoNSENsX result sheet</title>
  <meta name="description" content="COmpliance of NMR Structural
                                    ENSembles with eXperimental data" />
  <link rel="icon" href="http://consensx.itk.ppke.hu/public/_include/favicon.ico" type="image/x-icon">
  <link href="http://consensx.itk.ppke.hu/public/_include/css/bootstrap.min.css" rel="stylesheet">
  <link href="http://consensx.itk.ppke.hu/public/_include/css/style.css" rel="stylesheet">
  <link href="http://consensx.itk.ppke.hu/public/_include/css/main.css" rel="stylesheet">
  <link href='http://fonts.googleapis.com/css?family=Titillium+Web:400,200,200italic,300,300italic,400italic,600,600italic,700,700italic,900' rel='stylesheet' type='text/css'>
  <script src="http://consensx.itk.ppke.hu/public/_include/js/jquery.js"></script>
  <script src="http://consensx.itk.ppke.hu/public/_include/js/bootstrap.min.js"></script>
  <script src="http://consensx.itk.ppke.hu/public/_include/js/involve.js"></script>

  <!-- Add fancyBox -->
  <link rel="stylesheet" href="http://consensx.itk.ppke.hu/public/_include/fancybox/source/jquery.fancybox.css?v=2.1.5" type="text/css" media="screen" />
  <script type="text/javascript" src="http://consensx.itk.ppke.hu/public/_include/fancybox/source/jquery.fancybox.pack.js?v=2.1.5"></script>
  <link rel="stylesheet" href="http://consensx.itk.ppke.hu/public/_include/fancybox/source/helpers/jquery.fancybox-thumbs.css?v=1.0.7" type="text/css" media="screen" />
  <script type="text/javascript" src="http://consensx.itk.ppke.hu/public/_include/fancybox/source/helpers/jquery.fancybox-thumbs.js?v=1.0.7"></script>

  <script type="text/javascript">
    $(document).ready(function() {
      $(".fancybox").fancybox({
        openEffect  : 'elastic',
        closeEffect : 'elastic',
        openSpeed : 150,
        closeSpeed : 150,
        scrollOutside : true,
        });
    });
  </script>
</head>
<!--  PHP_VARIABLES -->""")

    html.close()


def writeFileTable(path, args, my_PDB, my_id, PDP_model_num,
                   my_NOE="not present", NOE_restraint_count=0):
    html = path + "result_sheet.php"
    html = open(html, 'a')

    if len(my_PDB.split('/')) != 1:
        my_PDB = my_PDB.split('/')[-1]

    if len(my_NOE.split('/')) != 1:
        my_NOE = my_NOE.split('/')[-1]

    if args.STR_file:
        if len(args.STR_file.split('/')) != 1:
            my_STR = args.STR_file.split('/')[-1]
        else:
            my_STR = args.STR_file
    else:
        my_STR = ""

    html.write("""
<body>
<div id="shortcodes" class="page" style="padding: 30px 0 0 0;">
    <div class="container" style="width: 1240px;">
    <div class="row" id="res-page-row" style="margin-left: 0px;">
      <div class="span12">
        <div class="title-page">
<h2 class="title"><b class="red"><a href="javascript:javascript:history.go(-1)"> CoNSEnsX<sup>+</sup></a></b></h2>
                    <h3 class="title-description"><b class="red">Co</b>mpliance of <b class="red">N</b>MR-derived <b class="red">S</b>tructural <b class="red">Ens</b>embles with e<b class="red">x</b>perimental data<br>+ selection</h3><br>
        <h3 class="title-description">Results sheet ID: <b class="red">{0}</b></h3>
        </div>
      </div>
    </div>

    <table class="files_table">
      <tr>
        <td class="head-td">PDB file:</td>
        <td><i>{1}</i></td>
        <td>{2} models found</td>
      </tr>
      <tr>
        <td class="head-td">NOE restraint file:</td>
        <td><i>{3}</i></td>
        <td>{4} distance restraints found</td>
      </tr>
      <tr>
        <td class="head-td">BMRB file:</td>
        <td><i>{5}</i></td>
        <td></td>
      </tr>
    </table>
    <div class="container">

      <div style="margin-left: auto; margin-right: auto; width: 520px;">

      <a href="values.csv" download>
      <button class="btn btn-result" id="csv_download" type="button">
      Download calculated data (CSV)
      </button>
      </a>

      <button class="btn btn-result" id="selection" type="button" data-toggle="collapse" data-target="#collapseExample" aria-expanded="false" aria-controls="collapseExample">
      Toggle selection options
      </button>

      </div>

    <div class="collapse" id="collapseExample">
      <table class="sel_controls" border="0">
        <tr class="form_param"><td>Compliance Measure: &nbsp;</td>
        <td>
        <div class="btn-group group1" data-toggle="buttons">

          <label class="btn btn-primary mes">
            <input type="radio" name="measures" class="options" value="correlation" > correlation
          </label>

          <label class="btn btn-primary mes">
            <input type="radio" name="measures" class="options" value="q-value" > Q-value
          </label>

          <label class="btn btn-primary mes">
            <input type="radio" name="measures" class="options" value="rmsd" > RMSD
          </label>

        </div>
        </td>
        <td>Min. ensemble size: &nbsp;</td>
            <td><input type="text" id="min_ens_size" value="2" name="sel_min_size"></td>

        </tr>
        <tr class="form_param"><td>Overdrive: &nbsp;</td>
            <td><input type="text" id="overdrive" name="sel_overdrive" value="5"></td>
            <td>Max. ensemble size: &nbsp;</td>
            <td><input type="text" id="max_ens_size" name="sel_max_size" value="{2}"></td>
        </tr>

      </table>
        <form name="sel_form" action="" method="post">
          <input type="text" id="sel_str" name="sel_string" style="visibility: hidden;">

          <button class="btn btn-primary" id="startsel" type="submit">
            Start selection
          </button>
        </form>

    </div>
  </div>

<?php

    if($_SERVER['REQUEST_METHOD'] == 'POST') {{
        if (isset($_POST['sel_string'])) {{

            $file = 'user_selection.txt';
            $command = $_POST['sel_string'];
            $command = str_replace("/", "\\n", $command);

            file_put_contents($file, $command);

            $my_dir = getcwd();
            $output = exec('../../consensx/selection.py ' . $my_dir . ' 2>&1');

            $sel_count = sizeof(explode(" ", $output));
            $sel_items = array();
            $measure   = "";
            $lines     = file($file);

            foreach ($lines as $line_num => $line) {{
                $tokens = explode(" ", $line);

                if ($tokens[0] == "MEASURE") {{
                    $measure = trim($tokens[1]);
                    continue;
                }} elseif ($tokens[0] == "OVERDRIVE") {{
                    continue;
                }} elseif ($tokens[0] == "MAX_SIZE") {{
                    continue;
                }} elseif ($tokens[0] == "MIN_SIZE") {{
                    continue;
                }} else {{
                    array_push($sel_items, $tokens[0]);
                }}

            }}

            $selection_data = 'selection_calced.txt';
            $lines = file($selection_data);
            $sel_data_array = array();

            foreach ($lines as $line_num => $line) {{
                $tokens = explode(" ", $line);
                $sel_data_array[$tokens[0]] = $tokens[1];
            }}
        }}

        echo '<h3 style="text-align: center;">Selection results</h3>';
        echo '<table class="files_table">';
        echo '<tr>';
        echo '<td class="head-td">' . strtoupper($measure) . '</td>';
        echo '<td>All models: ' . $model_count . '</td>';
        echo '<td>Selected models: ' . $sel_count .' <a href="selected.pdb" download>(download)</a></td>';
        echo '</tr>';

        foreach ($sel_items as $sel_item) {{
            echo '<tr>';
            echo '<td class="head-td" style="text-align: center;">' . $sel_item . '</td>';
            echo '<td style="text-align: center;">' . $full_pop[$measure][$sel_item] . '</td>';
            echo '<td style="text-align: center;">' . $sel_data_array[$sel_item] . '</td>';
            echo '</tr>';

        }}

        echo '</table>';
        echo '<h3 style="text-align: center;">Back-calculated data for the <b>original</b> ensemble</h3>';
    }}

?>""".format(my_id, my_PDB, PDP_model_num,
                   my_NOE, NOE_restraint_count,
                   my_STR))

    html.close()


def writeRDC_table_open(path, name, RDC_list_num):
    html = path + "result_sheet.php"
    html = open(html, 'a')

    html.write("""
    <div class="results">
      <h4 class="table-source">{0} {1}</h4>""".format(name, RDC_list_num))
    html.close()


def writeRDC_table_close(path):
    html = path + "result_sheet.php"
    html = open(html, 'a')

    html.write("    </div>\n")
    html.close()


def writeRDC_data(path, RDC_type, used_values, correl, q_value, rmsd,
                  corr_graph_name, graph_name, mod_corr_graph_name,
                  input_id):
    html = path + "result_sheet.php"
    html = open(html, 'a')

    html.write("""
      <table class="result_table">
        <tbody>
          <tr>
          <td><table class="values_table">
          <tr><td><strong>{0}</strong></td><td></td></tr>
          <tr><td>Values:</td><td>{1}</td></tr>
          <tr><td>Correlation:</td><td>{2}</td></tr>
          <tr><td>Q-factor:</td><td>{3} %</td></tr>
          <tr><td>RMSD:</td><td>{4}</td></tr>
          </table></td>
          <td><a class="fancybox" title="coming soon" href="{5}"><img width="270" src="{5}"></a></td>
          <td><a class="fancybox" title="coming soon" href="{6}"><img width="450" src="{6}"></a></td>
          <td><a class="fancybox" title="coming soon" href="{7}"><img width="270" src="{7}"></a></td>
          </tr>
        </tbody>
      </table>\n
      <div class="involve involve-result">
      <form oninput="x.value=parseInt({8}.value)">
        <table class="selection-table">
          <tr>
            <td class="involve-table"><label>Include with weight:</label></td>
            <td class="involve-table" style="width: 48px;"><output name="x" for="{8}">0</output></td>
            <td class="involve-table"><input type="range" class="inputrange" id="{8}" value="0" min="0" max="10"></td>
          </tr>
        </table>
        <p class="SVDnote hidden"><strong>SVD note:</strong> all RDCs of a group will have the same weight!</p>
        <p class="QVALnote hidden"><strong>Q-value as comp. measure:</strong> only RDCs are enabled!</p>
      </form>
      </div>\n""".format(RDC_type, used_values,
                            '{0:.3f}'.format(correl),
                            '{0:.3f}'.format(q_value),
                            '{0:.3f}'.format(rmsd),
                             corr_graph_name, graph_name, mod_corr_graph_name,
                             input_id))

    html.close()


def writeCS_data(path, RDC_type, used_values, correl, q_value, rmsd,
                  corr_graph_name, graph_name, mod_corr_graph_name,
                  input_id):
    html = path + "result_sheet.php"
    html = open(html, 'a')

    html.write("""
      <table class="result_table">
        <tbody>
          <tr>
          <td><table class="values_table">
          <tr><td><strong>{0}</strong></td><td></td></tr>
          <tr><td>Values:</td><td>{1}</td></tr>
          <tr><td>Correlation:</td><td>{2}</td></tr>
          <!-- <tr><td>Q-factor:</td><td>{3} %</td></tr> -->
          <tr><td>RMSD:</td><td>{4}</td></tr>
          </table></td>
          <td><a class="fancybox" title="coming soon" href="{5}"><img width="270" src="{5}"></a></td>
          <td><a class="fancybox" title="coming soon" href="{6}"><img width="450" src="{6}"></a></td>
          <td><a class="fancybox" title="coming soon" href="{7}"><img width="270" src="{7}"></a></td>
          </tr>
        </tbody>
      </table>\n
      <div class="involve involve-result">
      <form oninput="x.value=parseInt({8}.value)">
        <table class="selection-table">
          <tr>
            <td class="involve-table"><label>Include with weight:</label></td>
            <td class="involve-table" style="width: 48px;"><output name="x" for="{8}">0</output></td>
            <td class="involve-table"><input type="range" class="inputrange" id="{8}" value="0" min="0" max="10"></td>
          </tr>
        </table>
      </form>
      </div>\n""".format(RDC_type, used_values,
                            '{0:.3f}'.format(correl),
                            '{0:.3f}'.format(q_value),
                            '{0:.3f}'.format(rmsd),
                             corr_graph_name, graph_name, mod_corr_graph_name,
                             input_id))

    html.close()


def write_table_open(path, table_name):
    html = path + "result_sheet.php"
    html = open(html, 'a')

    html.write("""
    <div class="results">
      <h4 class="table-source">{0}</h4>""".format(table_name))
    html.close()


def write_table_close(path):
    html = path + "result_sheet.php"
    html = open(html, 'a')

    html.write("    </div>\n")
    html.close()


def write_table_data(path, data_type, used_values, correl, q_value, rmsd,
                     corr_graph_name, graph_name, input_id,
                     mod_corr_graph_name=None):
    html = path + "result_sheet.php"
    html = open(html, 'a')

    if mod_corr_graph_name:
      html.write("""
        <table class="result_table">
          <tbody>
            <tr>
            <td><table class="values_table">
            <tr><td><strong>{0}</strong></td><td></td></tr>
            <tr><td>Values:</td><td>{1}</td></tr>
            <tr><td>Correlation:</td><td>{2}</td></tr>
            <!-- <tr><td>Q-factor:</td><td>{3} %</td></tr> -->
            <tr><td>RMSD:</td><td>{4}</td></tr>
            </table></td>
            <td><a class="fancybox" title="coming soon" href="{5}"><img width="270" src="{5}"></a></td>
            <td><a class="fancybox" title="coming soon" href="{6}"><img width="450" src="{6}"></a></td>
            <td><a class="fancybox" title="coming soon" href="{7}"><img width="270" src="{7}"></a></td>
            </tr>
          </tbody>
        </table>\n
        <div class="involve involve-result">
        <form oninput="x.value=parseInt({8}.value)">
          <table class="selection-table">
            <tr>
              <td class="involve-table"><label>Include with weight:</label></td>
              <td class="involve-table" style="width: 48px;"><output name="x" for="{8}">0</output></td>
              <td class="involve-table"><input type="range" class="inputrange" id="{8}" value="0" min="0" max="10"></td>
            </tr>
          </table>
        </form>
        </div>\n""".format(data_type, used_values,
                              '{0:.3f}'.format(correl),
                              '{0:.3f}'.format(q_value),
                              '{0:.3f}'.format(rmsd),
                               corr_graph_name, graph_name,
                               mod_corr_graph_name, input_id))
    else:
        html.write("""
        <table class="result_table">
          <tbody>
            <tr>
            <td><table class="values_table">
            <tr><td><strong>{0}</strong></td><td></td></tr>
            <tr><td>Values:</td><td>{1}</td></tr>
            <tr><td>Correlation:</td><td>{2}</td></tr>
            <!-- <tr><td>Q-factor:</td><td>{3} %</td></tr> -->
            <tr><td>RMSD:</td><td>{4}</td></tr>
            </table></td>
              <td><a class="fancybox" title="coming soon" href="{5}"><img width="270" src="{5}"></a></td>
              <td><a class="fancybox" title="coming soon" href="{6}"><img width="450" src="{6}"></a></td>
            </tr>
          </tbody>
        </table>\n
        <div class="involve involve-result">
        <form oninput="x.value=parseInt({7}.value)">
          <table class="selection-table">
            <tr>
              <td class="involve-table"><label>Include with weight:</label></td>
              <td class="involve-table" style="width: 48px;"><output name="x" for="{7}">0</output></td>
              <td class="involve-table"><input type="range" class="inputrange" id="{7}" value="0" min="0" max="10"></td>
            </tr>
          </table>
        </form>
        </div>\n""".format(data_type, used_values,
                              '{0:.3f}'.format(correl),
                              '{0:.3f}'.format(q_value),
                              '{0:.3f}'.format(rmsd),
                               corr_graph_name, graph_name, input_id))

    html.close()


def write_sidechain_data(path, data_type, used_values, correl, q_value, rmsd,
                         input_id):
    html = path + "result_sheet.php"
    html = open(html, 'a')

    html.write("""
      <table class="result_table">
        <tbody>
          <tr>
          <td><table class="values_table">
          <tr><td><strong>{0}</strong></td><td></td></tr>
          <tr><td>Values:</td><td>{1}</td></tr>
          <tr><td>Correlation:</td><td>{2}</td></tr>
          <!-- <tr><td>Q-factor:</td><td>{3} %</td></tr> -->
          <tr><td>RMSD:</td><td>{4}</td></tr>
          </table></td>
          </tr>
        </tbody>
      </table>\n
      </div>\n""".format(data_type, used_values,
                            '{0:.3f}'.format(correl),
                            '{0:.3f}'.format(q_value),
                            '{0:.3f}'.format(rmsd),
                            input_id))


def write_bottom_table(path, NOE_violations, PRIDE_data):
    html = path + "result_sheet.php"
    html = open(html, 'a')

    html.write("""
    <div class="results">
      <h4 class="table-source">NOE violations and PRIDE-NMR</h4>
      <table class="result_table">

          <td><table class="NOE_PRIME_table">
          <tr><td><strong>NOE distance violation</strong></td></tr>
          <tr><td>Total # of violations:</td><td>{0}</td></tr>
          </table></td>
          <td width="10"></td>
          <td><a class="fancybox" title="NOE violations" href="NOE_hist.svg"><img width="320" src="NOE_hist.svg"></a></td>

          <td width="30"></td>

          <td><table class="NOE_PRIME_table">
          <tr><td><strong>PRIDE-NMR</strong></td></tr>
          <tr><td>Model with best score:</td><td>{1}</td></tr>
          <tr><td>Model with worst score:</td><td>{2}</td></tr>
          <tr><td>Average score:</td><td>{3}</td></tr>
          <tr><td>Standard Deviation:</td><td>{4}</td></tr>
          </table></td>
          <td width="10"></td>
          <td><a class="fancybox" title="NOE violations" href="PRIDE-NMR_score.svg"><img width="320" src="PRIDE-NMR_score.svg"></a></td>

      </table>
    </div>\n""".format(NOE_violations,
                       PRIDE_data[0],
                       PRIDE_data[1],
                       '{0:.3f}'.format(PRIDE_data[2]),
                       '{0:.3f}'.format(PRIDE_data[3])))

    html.close()


def close_HTML(path):
    html = path + "result_sheet.php"
    html = open(html, 'a')

    html.write("""
  </div>
</div>
</body>
</html>""")

    html.close()


def add_PHP_variables(path, PHP_dict, model_count):

    model_count = str(model_count)

    var_string = "<?php\n$model_count = " + model_count + ";\n"

    var_string += "$full_pop = array(\n"

    for measure in ["corr", "qval", "rmsd"]:
        if measure == "corr":
            array_name = "correlation"
        elif measure == "qval":
            array_name = "q-value"
        else:
            array_name = "rmsd"

        var_string += "'" + array_name + "' =>  array(\n"

        for key in PHP_dict.keys():
            if key.split('_')[-1] == measure:
                name = '_'.join(key.split('_')[:-1])
                var_string += "'" + name + "' => '" + PHP_dict[key] + "',\n"

        var_string += "),\n"

    var_string += ");\n?>\n"


    html = path + "result_sheet.php"

    f = open(html,'r')
    filedata = f.read()
    f.close()

    newdata = filedata.replace("<!--  PHP_VARIABLES -->", var_string)

    f = open(html,'w')
    f.write(newdata)
    f.close()

