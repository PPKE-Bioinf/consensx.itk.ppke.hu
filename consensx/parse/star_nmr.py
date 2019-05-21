import consensx.nmrpystar as nmrpystar
import pynmrstar


class RdcRecord:
    """Class for storing RDC data"""
    def __init__(self, resnum1, atom1, resnum2, atom2, rdc_value):
        self.RDC_type = (
            str(int(resnum1) - int(resnum2)) + '_' + atom1 + '_' + atom2
        )
        self.resnum = int(resnum1)
        self.atom = atom1
        self.resnum2 = int(resnum2)
        self.atom2 = atom2
        self.value = float(rdc_value)


class S2Record:
    """Class for storing S2 data"""
    def __init__(self, resnum, s2_type, s2_value):
        self.resnum = int(resnum)
        self.type = s2_type
        self.value = float(s2_value)
        self.calced = None


class JCoupRecord:
    """Class for storing J-Coupling data"""
    def __init__(self, resnum, jcoup_type, jcoup_value):
        self.resnum = int(resnum)
        self.type = jcoup_type
        self.value = float(jcoup_value)


class ChemShiftRecord:
    """Class for storing chemical shift data"""
    def __init__(self, resnum, res_name, atom_name, chemshift_value):
        self.resnum = int(resnum)
        self.res_name = res_name
        self.atom_name = atom_name
        self.value = float(chemshift_value)

    def __str__(self):
        return (
            "resnum: " + str(self.resnum) + "," +
            "res_name: " + self.res_name + "," +
            "atom_name: " + self.atom_name + "," +
            "value: " + str(self.value) + "\n"
        )

    def __repr__(self):
        return self.__str__()


class StarNMR:
    def __init__(self, str_file):
        self.parsed = pynmrstar.Entry.from_file(str_file)

    def parse_rdc(self):
        """Returns RDC lists as dictionaries containing RDC_Record objects,
        grouped by RDC type (keys())"""

        tag_list = ["Seq_ID_1", "Atom_ID_1", "Seq_ID_2", "Atom_ID_2", "Val"]
        rdc_lists = []
        rdc_records = []

        for rdc_loop in self.parsed.get_loops_by_category("RDC"):
            for row_data in rdc_loop.get_tag(tag_list):
                # ['95', 'H', '95', 'N', '-0.823']
                rdc_records.append(RdcRecord(*row_data))

            rdc_lists.append(rdc_records)

        # split list into dict according to RDC types
        new_rdc_list = []
        for RDC_list in rdc_lists:
            prev_type = ""
            rdc_dict = {}

            for record in RDC_list:
                if prev_type != record.RDC_type:
                    rdc_dict[record.RDC_type] = []
                    rdc_dict[record.RDC_type].append(record)
                else:
                    rdc_dict[record.RDC_type].append(record)

                prev_type = record.RDC_type

            new_rdc_list.append(rdc_dict)

        return new_rdc_list

    def parse_s2(self):
        """Returns a dictionary with the parsed S2 data"""

        try:
            saveShifts = self.parsed.value.saves["order_param"]

            loopShifts = saveShifts.loops[-1]
            s2_records = []

            for ix in range(len(loopShifts.rows)):  # fetch values from file
                row = loopShifts.getRowAsDict(ix)

                S2_value = float("{0:.2f}".format(float(row["S2_value"])))

                s2_records.append(
                    S2Record(
                        row["Residue_seq_code"], row["Atom_name"], S2_value)
                )

            # split list into dict according to S2 types
            S2_dict = {}
            prev_type = ""

            for record in s2_records:
                if prev_type != record.type:
                    S2_dict[record.type] = []
                    S2_dict[record.type].append(record)
                else:
                    S2_dict[record.type].append(record)

                prev_type = record.type

            return S2_dict

        except KeyError:
            print("No S2 parameter list found")
            return None

    def parse_s2_sidechain(self):
        """Returns a dictionary with the parsed S2 data"""

        sidechain_name = "side-chain_methyl_order_parameters"
        try:
            saveShifts = self.parsed.value.saves[sidechain_name]

            loopShifts = saveShifts.loops[-1]
            s2_records = []

            for ix in range(len(loopShifts.rows)):   # fetch values from file
                row = loopShifts.getRowAsDict(ix)

                S2_value = float("{0:.2f}".format(float(row["S2_value"])))

                s2_records.append(
                    S2Record(
                        row["Residue_seq_code"], row["Atom_name"], S2_value
                    )
                )

            return s2_records

        except KeyError:
            return None

    def parse_jcoup(self):
        """Returns a dictionary with the parsed J-coupling data"""
        try:
            saveShifts = self.parsed.value.saves["coupling_constants"]

            loopShifts = saveShifts.loops[-1]
            jcoup_records = []

            for ix in range(len(loopShifts.rows)):   # fetch values from file
                row = loopShifts.getRowAsDict(ix)

                JC_value = row["Coupling_constant_value"]
                JC_value = float("{0:.2f}".format(float(JC_value)))

                jcoup_records.append(
                    JCoupRecord(
                        row["Atom_one_residue_seq_code"],
                        row["Coupling_constant_code"],
                        JC_value
                    )
                )

            # split list into dict according to J-cuopling types
            jcoup_dict = {}
            prev_type = ""

            for record in jcoup_records:
                if prev_type != record.type:
                    jcoup_dict[record.type] = []
                    jcoup_dict[record.type].append(record)
                else:
                    jcoup_dict[record.type].append(record)

                prev_type = record.type

            return jcoup_dict

        except KeyError:
            print("No J-coupling parameter list found")
            return None

    def parse_chemshift(self):
        """Returns ChemShift lists as dictionaries containing ChemShift_Record
        objects, grouped by Atom_name (keys())"""

        chemshift_types = ["HA", "CA", "CB", "N", "H", "C"]
        tag_list = ["Seq_ID", "Comp_ID", "Atom_ID", "Atom_type", "Val"]
        chemshift_records = []
        chem_shift_lists = []
        ha_sum = 0.0

        for cs_loop in self.parsed.get_loops_by_category("Atom_chem_shift"):
            for row_data in cs_loop.get_tag(tag_list):
                row = dict(zip(tag_list, row_data))
                # ['21', 'PRO', 'HD3', 'H', '3.448']

                cs_value = float("{0:.2f}".format(float(row["Val"])))

                if row["Atom_ID"] == "HA2":
                    ha_sum += cs_value
                    continue

                if row["Atom_ID"] == "HA3":
                    ha_sum += cs_value

                    chemshift_records.append(
                        ChemShiftRecord(
                            row["Seq_ID"],
                            row["Comp_ID"],
                            "HA", ha_sum / 2
                        )
                    )
                    ha_sum = 0.0
                    continue

                if row["Atom_ID"] in chemshift_types:
                    chemshift_records.append(
                        ChemShiftRecord(
                            row["Seq_ID"],
                            row["Comp_ID"],
                            row["Atom_ID"],
                            cs_value
                        )
                    )

                elif row["Atom_ID"] == "HN":
                    chemshift_records.append(
                        ChemShiftRecord(
                            row["Seq_ID"],
                            row["Comp_ID"],
                            "H",
                            cs_value
                        )
                    )

            chem_shift_lists.append(chemshift_records)

        new_cs_list = []

        for chem_shift_list in chem_shift_lists:
            chem_shift_dict = {}

            for record in chem_shift_list:

                if record.atom_name in list(chem_shift_dict.keys()):
                    chem_shift_dict[record.atom_name].append(record)
                else:
                    chem_shift_dict[record.atom_name] = []
                    chem_shift_dict[record.atom_name].append(record)

            new_cs_list.append(chem_shift_dict)

        return new_cs_list
