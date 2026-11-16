"""Object-scoped resname discovery, filtering and selection shortcuts."""
import os
os.environ.setdefault('QT_QPA_PLATFORM','offscreen')
import unittest
from pymol import cmd
from pymol.Qt import QtCore, QtWidgets
from pmg_qt.visual_representations_gui import VisualRepresentationsDialog


class ResnameShortcutTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app=QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def setUp(self):
        cmd.reinitialize()
        for name,resnames in [('first_source',['AL31','AL31','DPPC','U']),
                              ('second_source',['CHOL','CHOL','POPC'])]:
            for i,resn in enumerate(resnames):
                cmd.pseudoatom(name,name='C'+str(i),resn=resn,resi=i+1,pos=[i,0,0])
        self.parent=QtWidgets.QWidget();self.parent.cmd=cmd
        self.dialog=VisualRepresentationsDialog(self.parent,cmd)
        self.dialog.select_target('first_source')

    def tearDown(self):
        self.dialog.close();self.dialog.deleteLater();self.parent.deleteLater()
        self.app.processEvents()

    def names(self):
        listing=self.dialog.list_resnames
        return [listing.item(i).data(QtCore.Qt.ItemDataRole.UserRole) for i in range(listing.count())]

    def test_unique_names_click_count_and_scope(self):
        d=self.dialog
        self.assertEqual(self.names(),['AL31','DPPC','U'])
        item=d.list_resnames.item(0)
        d._choose_resname(item)
        self.assertEqual(d.ed_sele.text(),'resname AL31')
        self.assertIn('2',d.lbl_preview.text())
        self.assertFalse(d.ctrl.owned)
        self.assertIn('first_source',cmd.get_names('objects',enabled_only=1))
        self.assertIn('second_source',cmd.get_names('objects',enabled_only=1))

    def test_search_switch_and_refresh_replace_old_data(self):
        d=self.dialog;d.ed_resname_search.setText('dp')
        shown=[d.list_resnames.item(i).data(QtCore.Qt.ItemDataRole.UserRole)
               for i in range(d.list_resnames.count()) if not d.list_resnames.item(i).isHidden()]
        self.assertEqual(shown,['DPPC'])
        d.select_target('second_source')
        self.assertEqual(self.names(),['CHOL','POPC'])
        self.assertEqual(d.ed_resname_search.text(),'')
        cmd.pseudoatom('second_source',resn='NEW',pos=[8,0,0])
        d.btn_refresh.click()
        self.assertEqual(self.names(),['CHOL','NEW','POPC'])

    def test_auto_apply_and_manual_expression_are_preserved(self):
        d=self.dialog
        d.cmb_style.setCurrentIndex(d.cmb_style.findData('spheres'))
        d.chk_auto.setChecked(True)
        d._choose_resname(d.list_resnames.item(1))
        self.assertTrue(d.ctrl.last_report.ok,d.status.text())
        self.assertEqual(cmd.count_atoms(d.ctrl.owned[0]),1)
        self.assertIn('second_source',cmd.get_names('objects',enabled_only=1))
        d.chk_auto.setChecked(False)
        d.ed_sele.setText('resname DPPC or nucleic')
        d._on_sele_finished()
        self.assertEqual(d.ctrl.layers[0].selection,'resname DPPC or nucleic')
        d.restore_source()

    def test_many_names_scroll_without_growing_the_window(self):
        d=self.dialog
        for i in range(80):
            cmd.pseudoatom('second_source',resn='R%03d'%i,pos=[i,0,0])
        d.select_target('second_source');d.resize(760,860);d.show()
        self.app.processEvents()
        self.assertEqual(d.list_resnames.count(),82)
        self.assertLessEqual(d.list_resnames.height(),110)
        self.assertGreater(d.list_resnames.verticalScrollBar().maximum(),0)
        d.ed_resname_search.setText('R079')
        self.assertEqual(sum(not d.list_resnames.item(i).isHidden()
                             for i in range(d.list_resnames.count())),1)


if __name__=='__main__': unittest.main()
