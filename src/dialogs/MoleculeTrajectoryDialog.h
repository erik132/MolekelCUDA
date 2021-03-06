#ifndef MOLECULETRAJECTORYDIALOG_H_
#define MOLECULETRAJECTORYDIALOG_H_
//
// Molekel - Molecular Visualization Program
// Copyright (C) 2006, 2007, 2008, 2009 Swiss National Supercomputing Centre (CSCS)
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA  02110-1301, USA.
//
// $Author$
// $Date$
// $Revision$
//

// QT
#include <QWidget>
#include <QDialog>
#include <QVBoxLayout>
#include <QPushButton>


#include "../widgets/MoleculeTrajectoryWidget.h"

class MolekelMolecule;
class AtomTrajectoryAnimator;

//
//-----------------------
//   -----------------
//  |                |
//  |                |
//  |                |
//  |    widget      |
//  |                |
//  |                |
//  |                |
//   ----------------
//   ____
//  |_Ok_|
//-----------------------
//

/// Container for MoleculeTrajectoryWidget
class MoleculeTrajectoryDialog : public QDialog
{
    Q_OBJECT
public:
    /// Constructor.
    /// @param mol reference to molecule to edit
    /// @param parent widget
    MoleculeTrajectoryDialog( MolekelMolecule *mol,
                              AtomTrajectoryAnimator* a,
                              QWidget* parent = 0 ) : QDialog( parent )
    {
        QVBoxLayout* mainLayout = new QVBoxLayout();
        MoleculeTrajectoryWidget* mvw = new MoleculeTrajectoryWidget( this );
        mvw->SetMolecule( mol );
        mvw->SetAnimator( a );
        mainLayout->addWidget( mvw );
        QPushButton* b = new QPushButton( tr( "Close" ) );
        connect( b, SIGNAL( released() ), this, SLOT( AcceptSlot() ) );
        mainLayout->addWidget( b );
        this->setLayout( mainLayout );
    }
public slots:
    /// Called when 'Ok' clicked.
    void AcceptSlot() { accept(); }
};

#endif /*MOLECULETRAJECTORYDIALOG_H_*/
