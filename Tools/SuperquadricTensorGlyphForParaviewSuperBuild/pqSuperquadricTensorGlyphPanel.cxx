#include "pqSuperquadricTensorGlyphPanel.h"

pqSuperquadricTensorGlyphPanel::pqSuperquadricTensorGlyphPanel(pqProxy *proxy, QWidget *p)
  : Superclass(proxy, p)
{
  // -------------------------------------------------------------------------
  // Automatically disable gamma roundness and enable both remaining
  // roundness widgets if roundness should not be set individually per tensor.

  this->FixedThetaPhiRoundnessWidget = this->findChild<QCheckBox*>("FixedThetaPhiRoundness");
  if (! this->FixedThetaPhiRoundnessWidget)
    {
    qWarning() << "Failed to locate FixedThetaPhiRoundness widget.";
    return;
    }

  this->GammaRoundnessWidget = this->findChild<QWidget*>("GammaRoundness");
  if (this->GammaRoundnessWidget)
    QObject::connect(this->FixedThetaPhiRoundnessWidget, SIGNAL(toggled(bool)),
                     this->GammaRoundnessWidget, SLOT(setDisabled(bool)));

  this->GammaRoundnessLabel = this->findChild<QLabel*>("_labelForGammaRoundness");
  if (this->GammaRoundnessLabel)
    QObject::connect(this->FixedThetaPhiRoundnessWidget, SIGNAL(toggled(bool)),
                     this->GammaRoundnessLabel, SLOT(setDisabled(bool)));


  this->ThetaRoundnessWidget = this->findChild<QWidget*>("ThetaRoundness");
  if (this->ThetaRoundnessWidget)
    QObject::connect(this->FixedThetaPhiRoundnessWidget, SIGNAL(toggled(bool)),
                     this->ThetaRoundnessWidget, SLOT(setEnabled(bool)));

  this->ThetaRoundnessLabel = this->findChild<QLabel*>("_labelForThetaRoundness");
  if (this->ThetaRoundnessLabel)
    QObject::connect(this->FixedThetaPhiRoundnessWidget, SIGNAL(toggled(bool)),
                     this->ThetaRoundnessLabel, SLOT(setEnabled(bool)));


  this->PhiRoundnessWidget = this->findChild<QWidget*>("PhiRoundness");
  if (this->PhiRoundnessWidget)
    QObject::connect(this->FixedThetaPhiRoundnessWidget, SIGNAL(toggled(bool)),
                     this->PhiRoundnessWidget, SLOT(setEnabled(bool)));

  this->PhiRoundnessLabel = this->findChild<QLabel*>("_labelForPhiRoundness");
  if (this->PhiRoundnessLabel)
    QObject::connect(this->FixedThetaPhiRoundnessWidget, SIGNAL(toggled(bool)),
                     this->PhiRoundnessLabel, SLOT(setEnabled(bool)));

  this->FixedThetaPhiRoundnessWidget->toggle();
  this->FixedThetaPhiRoundnessWidget->toggle();


  // -------------------------------------------------------------------------
  // Automatically disable "color by" dropdown box if ColorGlyphs is not set

  this->ColorGlyphsWidget = this->findChild<QCheckBox*>("ColorGlyphs");
  if (! this->ColorGlyphsWidget)
    {
    qWarning() << "Failed to locate ColorGlyphs widget.";
    return;
    }

  this->ColorModeWidget = this->findChild<QWidget*>("ColorMode");
  if (this->ColorModeWidget)
    QObject::connect(this->ColorGlyphsWidget, SIGNAL(toggled(bool)),
                     this->ColorModeWidget, SLOT(setEnabled(bool)));

  this->ColorModeLabel = this->findChild<QLabel*>("_labelForColorMode");
  if (this->ColorModeLabel)
    QObject::connect(this->ColorGlyphsWidget, SIGNAL(toggled(bool)),
                     this->ColorModeLabel, SLOT(setEnabled(bool)));

  this->ColorGlyphsWidget->toggle();
  this->ColorGlyphsWidget->toggle();


  // -------------------------------------------------------------------------
  // Automatically disable "MaxScaleFactor" textbox if LimitScaling is not set

  this->LimitScalingByEigenvaluesWidget = this->findChild<QCheckBox*>("LimitScalingByEigenvalues");
  if (! this->LimitScalingByEigenvaluesWidget)
    {
    qWarning() << "Failed to locate LimitScalingByEigenvalues widget.";
    return;
    }

  this->MaxScaleFactorWidget = this->findChild<QWidget*>("MaxScaleFactor");
  if (this->MaxScaleFactorWidget)
    QObject::connect(this->LimitScalingByEigenvaluesWidget, SIGNAL(toggled(bool)),
                     this->MaxScaleFactorWidget, SLOT(setEnabled(bool)));

  this->MaxScaleFactorLabel = this->findChild<QLabel*>("_labelForMaxScaleFactor");
  if (this->MaxScaleFactorLabel)
    QObject::connect(this->LimitScalingByEigenvaluesWidget, SIGNAL(toggled(bool)),
                     this->MaxScaleFactorLabel, SLOT(setEnabled(bool)));

  this->LimitScalingByEigenvaluesWidget->toggle();
  this->LimitScalingByEigenvaluesWidget->toggle();


  // -------------------------------------------------------------------------
  // Automatically switch to color by eigenvalues if input data set has no
  // scalars.
  QComboBox* combo1 = this->findChild<QComboBox*>("SelectInputScalars");
  if (combo1)
    {
    if (combo1->count() < 1)
      {
      QComboBox* combo2 = this->findChild<QComboBox*>("ColorMode");
      if (combo2)
        {
        // Locate index of item whose text is "eigenvalues"
        for (unsigned int i = 0; i < combo2->count(); i++)
          {
          if (combo2->itemText(i) == QString("eigenvalues"))
            {
            combo2->setCurrentIndex(i);
            break;
            }
          }
        }
      }
    }

  return;
}
