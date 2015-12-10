import javax.swing.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.table.DefaultTableModel;
import java.awt.event.*;
import java.util.HashSet;

public class dlgAddCategory extends JDialog {

    private JPanel contentPane;
    private JButton buttonOK;
    private JButton buttonCancel;
    private JTextField txtCategoryName;
    private JList lstCategoryValues;
    private JButton btnAdd;
    private JButton btnRemove;
    private JLabel lblCategoryName;
    private JTable tblCategoryValues;

    public dlgAddCategory() {
        setContentPane(contentPane);
        setModal(true);
        getRootPane().setDefaultButton(buttonOK);
        this.setTitle("Add experimental category");

        initCategoryModel();
        initEventListeners();

        // TODO Disable the OK button until at least a title and value have been added.

// call onCancel() when cross is clicked
        setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
        addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent e) {
                onCancel();
            }
        });

// call onCancel() on ESCAPE
        contentPane.registerKeyboardAction(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                onCancel();
            }
        }, KeyStroke.getKeyStroke(KeyEvent.VK_ESCAPE, 0), JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
    }

    // TODO Better name than not sure here.
    public enum DialogResult {
        NOTSURE,
        OK,
        CANCEL
    }

    public DialogResult Result = DialogResult.NOTSURE;

    public SampleCategory ProvidedCategory;

    private void initCategoryModel() {
        DefaultTableModel model = new DefaultTableModel();
        model.addColumn("Value");
        tblCategoryValues.setModel(model);
    }

    private void initEventListeners() {
        buttonOK.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                onOK();
            }
        });

        buttonCancel.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                onCancel();
            }
        });

        btnAdd.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                // TODO Ideally the row would be highlighted and the user could start editing immediately after pressing "Add".
                ((DefaultTableModel)tblCategoryValues.getModel()).addRow(new Object[] {"TODO Value here"});

                ValidateFields();
            }
        });

        // TODO Fix form scaling positioning.
        btnRemove.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                int[] categoryValues = tblCategoryValues.getSelectedRows();

                for (int categoryValue : categoryValues) {
                    ((DefaultTableModel)tblCategoryValues.getModel()).removeRow(categoryValue);
                }

                ValidateFields();
            }
        });

        txtCategoryName.getDocument().addDocumentListener(new DocumentListener() {
            @Override
            public void insertUpdate(DocumentEvent e) {
                ValidateFields();
            }

            @Override
            public void removeUpdate(DocumentEvent e) {
                ValidateFields();
            }

            @Override
            public void changedUpdate(DocumentEvent e) {

            }
        });
    }

    private void ValidateFields() {
        boolean isValid = true;

        if (txtCategoryName.getText().length() == 0) {
            // An experimental category without a name shouldn't be add-able.
            isValid = false;
        }

        if (tblCategoryValues.getRowCount() == 0) {
            // An experimental category without any defined values shouldn't be add-able.
            isValid = false;
        }

        buttonOK.setEnabled(isValid);
        if (isValid) {
            buttonOK.setToolTipText("");
        } else {
            // TODO Use constants here.
            // TODO Do something better than setting the tooltip text - it's not obvious enough.
            buttonOK.setToolTipText("Please provide a category title and at least one valid value.");
        }
    }

    private void onOK() {
        // TODO Just write a function for converting from the UI to an object as I'm doing here.
        HashSet<String> validValues = new HashSet<>();
        for (int value = 0; value < tblCategoryValues.getRowCount(); value++) {
            validValues.add((String) tblCategoryValues.getValueAt(value, 0));
        }

        ProvidedCategory = new SampleCategory(txtCategoryName.getText(),
                                              validValues);

        Result = DialogResult.OK;
        dispose();
    }

    private void onCancel() {
        Result = DialogResult.CANCEL;
        dispose();
    }

    public DialogResult Display() {
        pack();

        // TODO Make an extension method or similar to center a form.
        // TODO This is not working.
        // Center the dialog on-screen.
        setLocationRelativeTo(null);

        setVisible(true);

        return Result;
    }
}