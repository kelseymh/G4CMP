////////////////////////////////////////////////////////
//
// RISQTutorialDetectorParameters.hh
//
// linehan3@fnal.gov
// This is a file of detector parameters, those
// that define dimensions of the various geometries
// used in the silicon qubit sims code.
//
////////////////////////////////////////////////////////

#ifndef RISQTutorialDetectorParameters_h
#define RISQTutorialDetectorParameters_h 1

#include "CLHEP/Units/SystemOfUnits.h"
#include <cstdlib>
#include <cmath>

namespace RISQTutorialDetectorParameters
{
  //----------------------------------------------------------------
  //Overall World
  constexpr double dp_worldSize = 6.0 * CLHEP::cm;

  //Misc
  //----------------------------------------------------------------
  constexpr double dp_eps = 0.0001*CLHEP::mm;
  constexpr double pi = 3.141592654;
  
  //----------------------------------------------------------------
  //Silicon chip dimensions
  constexpr double dp_siliconChipDimX = 8.0 * CLHEP::mm;
  constexpr double dp_siliconChipDimY = 8.0 * CLHEP::mm;
  constexpr double dp_siliconChipDimZ = 0.381 * CLHEP::mm;


  //----------------------------------------------------------------
  //Parameters of the qubit housing (currently copper)
  constexpr bool dp_useQubitHousing = true;
  constexpr double dp_housingDimX = 16.0 * CLHEP::mm;
  constexpr double dp_housingDimY = 16.0 * CLHEP::mm;
  constexpr double dp_housingDimZ = 1.0 * CLHEP::cm; //Used to be 1
  constexpr double dp_housingCentralCutoutDimX = dp_siliconChipDimX + dp_eps;
  constexpr double dp_housingCentralCutoutDimY = dp_siliconChipDimY + dp_eps;
  constexpr double dp_housingCentralCutoutDimZ = dp_housingDimZ * 0.3;
  constexpr double dp_housingRadialCutoutDimX = 1.03 * CLHEP::mm; //Much of the cutout is actually thinner, but the important thing is the corner that the chip sits on (to leading order) - the wall-chip interfaces will be a smidge shorter than in reality, but for now this is fine.
  constexpr double dp_housingRadialCutoutDimY = 2 * CLHEP::mm; //Arbitrary, but 2 mm should do the trick
  constexpr double dp_housingRadialCutoutDimZ = dp_siliconChipDimZ; //QUBIT IS FLUSH WITH TOP OF HOUSING

  //----------------------------------------------------------------
  //Parameters of the niobium film ground plane
  constexpr double dp_useGroundPlane = true;
  constexpr double dp_groundPlaneDimX = dp_siliconChipDimX;
  constexpr double dp_groundPlaneDimY = dp_siliconChipDimY;
  constexpr double dp_groundPlaneDimZ = 90 * CLHEP::nm;

  //----------------------------------------------------------------
  //Parameters of the transmission line
  constexpr bool dp_useTransmissionLine = true;
  constexpr double dp_transmissionLineBaseLayerDimY = 626 * CLHEP::um;
  constexpr double dp_transmissionLineBaseLayerDimX = dp_groundPlaneDimX;
  constexpr double dp_transmissionLineBaseLayerDimZ = dp_groundPlaneDimZ;
  constexpr double dp_transmissionLineCavityFullWidth = 22 * CLHEP::um;
  constexpr double dp_transmissionLineConductorWidth = 10 * CLHEP::um;
  constexpr double dp_transmissionLinePad1Offset = -3562.5 * CLHEP::um;
  constexpr double dp_transmissionLinePad2Offset = 3562.5 * CLHEP::um;

  
  //----------------------------------------------------------------
  //Generic pad parameters. (Pad as defined here is arrow pointing to right)
  constexpr double dp_padEmptyPart1DimX = 500 * CLHEP::um;
  constexpr double dp_padEmptyPart1DimY = 625 * CLHEP::um;
  constexpr double dp_padEmptyPart1DimZ = dp_groundPlaneDimZ;
  
  //This part of the pad has to get rotated, unfortunately, so the X, Y, and Z parameters
  //are pre-rotation for clarity of definition, rather than consistency with the rest of the
  //variables here.
  constexpr double dp_padEmptyPart2TrdZ = 233.432 * CLHEP::um;
  constexpr double dp_padEmptyPart2TrdX1 = dp_padEmptyPart1DimY;
  constexpr double dp_padEmptyPart2TrdX2 = dp_transmissionLineCavityFullWidth;
  constexpr double dp_padEmptyPart2TrdY1 = dp_padEmptyPart1DimZ;
  constexpr double dp_padEmptyPart2TrdY2 = dp_padEmptyPart1DimZ;

  //Non-empty parts of the pad
  constexpr double dp_padPart1DimX = 375 * CLHEP::um;
  constexpr double dp_padPart1DimY = 375 * CLHEP::um;
  constexpr double dp_padPart1DimZ = dp_groundPlaneDimZ;
  
  //This part of the pad has to get rotated, unfortunately, so the X, Y, and Z parameters
  //are pre-rotation for clarity of definition, rather than consistency with the rest of the
  //variables here.
  constexpr double dp_padPart2TrdZ = 233.432 * CLHEP::um;
  constexpr double dp_padPart2TrdX1 = dp_padPart1DimY;
  constexpr double dp_padPart2TrdX2 = dp_transmissionLineConductorWidth;
  constexpr double dp_padPart2TrdY1 = dp_padPart1DimZ;
  constexpr double dp_padPart2TrdY2 = dp_padPart1DimZ;
  constexpr double dp_padPart2InternalShiftX = (dp_padEmptyPart1DimX-dp_padPart1DimX)/2.0;

  //----------------------------------------------------------------
  //Resonator assembly parameters. Y extent of the mother volume goes from the transmission line cavity side to
  //the flux line cavity side.
  constexpr bool dp_useResonatorAssembly = true;
  constexpr double dp_resonatorAssemblyBaseNbDimX = 891.618 * CLHEP::um;
  constexpr double dp_resonatorAssemblyBaseNbDimY = 1925.311 * CLHEP::um;
  constexpr double dp_resonatorAssemblyBaseNbDimZ = dp_groundPlaneDimZ;
  constexpr double dp_resonatorAssemblyBaseNbEdgeBottomDimY = 5.0 * CLHEP::um; //Separation between the transmission line empty edge and the nearby resonator empty edge
  
  constexpr double dp_tlCouplingEmptyDimX = 447.518 * CLHEP::um;
  constexpr double dp_tlCouplingEmptyDimY = 22 * CLHEP::um;
  constexpr double dp_tlCouplingEmptyDimZ = dp_groundPlaneDimZ;
  constexpr double dp_tlCouplingConductorDimX = dp_tlCouplingEmptyDimX;
  constexpr double dp_tlCouplingConductorDimY = 10 * CLHEP::um;
  constexpr double dp_tlCouplingConductorDimZ = dp_groundPlaneDimZ;

  //Curve parameters
  constexpr double dp_resonatorAssemblyCurveSmallestRadius = 88.763 * CLHEP::um;
  constexpr double dp_resonatorAssemblyCurveCentralRadius = dp_resonatorAssemblyCurveSmallestRadius + dp_tlCouplingEmptyDimY / 2.0;
  constexpr double dp_curveEmptyDimZ = dp_groundPlaneDimZ;

  //Straight horizontal line 1
  constexpr double dp_shl1EmptyDimX = 52.642*CLHEP::um;
  constexpr double dp_shl1EmptyDimY = dp_tlCouplingEmptyDimY;
  constexpr double dp_shl1EmptyDimZ = dp_groundPlaneDimZ;
  constexpr double dp_shl1ConductorDimX = dp_shl1EmptyDimX;
  constexpr double dp_shl1ConductorDimY = dp_tlCouplingConductorDimY;
  constexpr double dp_shl1ConductorDimZ = dp_groundPlaneDimZ;

  //No extra params needed for half circle 1
  
  //Straight horizontal line 2
  constexpr double dp_shl2EmptyDimX = 287.796*CLHEP::um;
  constexpr double dp_shl2EmptyDimY = dp_tlCouplingEmptyDimY;
  constexpr double dp_shl2EmptyDimZ = dp_groundPlaneDimZ;
  constexpr double dp_shl2ConductorDimX = dp_shl2EmptyDimX;
  constexpr double dp_shl2ConductorDimY = dp_tlCouplingConductorDimY;
  constexpr double dp_shl2ConductorDimZ = dp_groundPlaneDimZ;

  //No extra params needed for half circle 2
  
  //Straight horizontal line 3
  constexpr double dp_shl3EmptyDimX = dp_shl2EmptyDimX;
  constexpr double dp_shl3EmptyDimY = dp_tlCouplingEmptyDimY;
  constexpr double dp_shl3EmptyDimZ = dp_groundPlaneDimZ;
  constexpr double dp_shl3ConductorDimX = dp_shl3EmptyDimX;
  constexpr double dp_shl3ConductorDimY = dp_tlCouplingConductorDimY;
  constexpr double dp_shl3ConductorDimZ = dp_groundPlaneDimZ;

  //No extra params needed for half circle 3
  
  //Straight horizontal line 4
  constexpr double dp_shl4EmptyDimX = dp_shl3EmptyDimX;
  constexpr double dp_shl4EmptyDimY = dp_tlCouplingEmptyDimY;
  constexpr double dp_shl4EmptyDimZ = dp_groundPlaneDimZ;
  constexpr double dp_shl4ConductorDimX = dp_shl4EmptyDimX;
  constexpr double dp_shl4ConductorDimY = dp_tlCouplingConductorDimY;
  constexpr double dp_shl4ConductorDimZ = dp_groundPlaneDimZ;

  //Straight horizontal line 5
  constexpr double dp_shl5EmptyDimX = dp_shl3EmptyDimX;
  constexpr double dp_shl5EmptyDimY = dp_tlCouplingEmptyDimY;
  constexpr double dp_shl5EmptyDimZ = dp_groundPlaneDimZ;
  constexpr double dp_shl5ConductorDimX = dp_shl5EmptyDimX;
  constexpr double dp_shl5ConductorDimY = dp_tlCouplingConductorDimY;
  constexpr double dp_shl5ConductorDimZ = dp_groundPlaneDimZ;

  //Straight horizontal line 6
  constexpr double dp_shl6EmptyDimX = dp_shl3EmptyDimX;
  constexpr double dp_shl6EmptyDimY = dp_tlCouplingEmptyDimY;
  constexpr double dp_shl6EmptyDimZ = dp_groundPlaneDimZ;
  constexpr double dp_shl6ConductorDimX = dp_shl6EmptyDimX;
  constexpr double dp_shl6ConductorDimY = dp_tlCouplingConductorDimY;
  constexpr double dp_shl6ConductorDimZ = dp_groundPlaneDimZ;

  //Straight horizontal line 7
  constexpr double dp_shl7EmptyDimX = 44.705*CLHEP::um;
  constexpr double dp_shl7EmptyDimY = dp_tlCouplingEmptyDimY;
  constexpr double dp_shl7EmptyDimZ = dp_groundPlaneDimZ;
  constexpr double dp_shl7ConductorDimX = dp_shl7EmptyDimX;
  constexpr double dp_shl7ConductorDimY = dp_tlCouplingConductorDimY;
  constexpr double dp_shl7ConductorDimZ = dp_groundPlaneDimZ;


  //Straight vertical line 1
  constexpr double dp_svl1EmptyDimX = dp_tlCouplingEmptyDimY;
  constexpr double dp_svl1EmptyDimY = 47.611*CLHEP::um;
  constexpr double dp_svl1EmptyDimZ = dp_groundPlaneDimZ;
  constexpr double dp_svl1ConductorDimX = dp_tlCouplingConductorDimY;
  constexpr double dp_svl1ConductorDimY = dp_svl1EmptyDimY;
  constexpr double dp_svl1ConductorDimZ = dp_groundPlaneDimZ;


  //Shunt coupler
  constexpr double dp_shuntCouplerHorizontalEmptyDimX = 194.614*CLHEP::um;
  constexpr double dp_shuntCouplerHorizontalEmptyDimY = 38.206*CLHEP::um;
  constexpr double dp_shuntCouplerHorizontalEmptyDimZ = dp_groundPlaneDimZ;
  constexpr double dp_shuntCouplerHorizontalConductorDimX = 187.451*CLHEP::um;
  constexpr double dp_shuntCouplerHorizontalConductorDimY = 31.64*CLHEP::um;
  constexpr double dp_shuntCouplerHorizontalConductorDimZ = dp_groundPlaneDimZ;
  constexpr double dp_shuntCouplerHorizontalConductorNubDimX = dp_tlCouplingConductorDimY;
  constexpr double dp_shuntCouplerHorizontalConductorNubDimY = 0.5*(dp_shuntCouplerHorizontalEmptyDimY-dp_shuntCouplerHorizontalConductorDimY);
  constexpr double dp_shuntCouplerHorizontalConductorNubDimZ = dp_groundPlaneDimZ;
  constexpr double dp_shuntCouplerLobeEmptyDimX = 58.187*CLHEP::um;
  constexpr double dp_shuntCouplerLobeEmptyDimY = 58.187*CLHEP::um;
  constexpr double dp_shuntCouplerLobeEmptyDimZ = dp_groundPlaneDimZ;
  constexpr double dp_shuntCouplerLobeConductorDimX = 50.000*CLHEP::um;
  constexpr double dp_shuntCouplerLobeConductorDimY = 55.957*CLHEP::um;
  constexpr double dp_shuntCouplerLobeConductorDimZ = dp_groundPlaneDimZ;


  //Shunt itself
  constexpr double dp_shuntCenterToBottomRightCornerOfBaseLayerDimX = 555 *CLHEP::um;//617.5 * CLHEP::um;
  constexpr double dp_shuntCenterToBottomRightCornerOfBaseLayerDimY = 1762.118 * CLHEP::um;
  constexpr double dp_shuntHorizontalBlockEmptyDimX = 304. * CLHEP::um;
  constexpr double dp_shuntHorizontalBlockEmptyDimY = 72. * CLHEP::um;
  constexpr double dp_shuntHorizontalBlockEmptyDimZ = dp_groundPlaneDimZ;
  constexpr double dp_shuntVertBlockEmptyDimX = dp_shuntHorizontalBlockEmptyDimY;
  constexpr double dp_shuntVertBlockEmptyDimY = dp_shuntHorizontalBlockEmptyDimX;
  constexpr double dp_shuntVertBlockEmptyDimZ = dp_groundPlaneDimZ;
  constexpr double dp_shuntHorizontalBlockConductorDimX = 280. * CLHEP::um;
  constexpr double dp_shuntHorizontalBlockConductorDimY = 24. * CLHEP::um;
  constexpr double dp_shuntHorizontalBlockConductorDimZ = dp_groundPlaneDimZ;
  constexpr double dp_shuntVertBlockConductorDimX = dp_shuntHorizontalBlockConductorDimY;
  constexpr double dp_shuntVertBlockConductorDimY = 256.*CLHEP::um;
  constexpr double dp_shuntVertBlockConductorDimZ = dp_groundPlaneDimZ;

  
  //Overall resonator assembly properties
  constexpr double dp_resonatorLateralSpacing = 1800*CLHEP::um + 56*CLHEP::um; //Adding arbitrary 64 um to spacing to get things to line up.
  constexpr double dp_centralResonatorOffsetX = -1*dp_resonatorAssemblyBaseNbDimX/2.0 + 617.5*CLHEP::um - 62*CLHEP::um; //Needs to be defined because the center of the resonator object is not the center of the square. Subtracting an additional arbitrary 68 um to the offset to get things to line up. Good enough.
  
  



  //----------------------------------------------------------------
  //Flux line parameters

  constexpr bool dp_useFluxLines = true;
  
  //Straight flux line info
  constexpr double dp_fluxLineBaseNbLayerDimX = 625 * CLHEP::um;
  constexpr double dp_fluxLineBaseNbLayerDimY = 1882.005 * CLHEP::um;
  constexpr double dp_fluxLineBaseNbLayerDimZ = dp_groundPlaneDimZ;

  constexpr double dp_fluxLinePadOffsetY = 0.5*dp_fluxLineBaseNbLayerDimY - 0.5*dp_padEmptyPart1DimX;

  constexpr double dp_fluxLineEmptyDimX = 22 * CLHEP::um;
  constexpr double dp_fluxLineEmptyDimY = 1140.498 * CLHEP::um;
  constexpr double dp_fluxLineEmptyDimZ = dp_groundPlaneDimZ;
  constexpr double dp_fluxLineLineOffsetY = 0.5*dp_fluxLineBaseNbLayerDimY - dp_padEmptyPart1DimX - dp_padEmptyPart2TrdZ - 0.5*dp_fluxLineEmptyDimY;
  constexpr double dp_fluxLineConductorDimX = 10.15* CLHEP::um;
  constexpr double dp_fluxLineConductorDimY = dp_fluxLineEmptyDimY;
  constexpr double dp_fluxLineConductorDimZ = dp_groundPlaneDimZ;


  //Corner flux line info (initially made for pad closest to top left corner of qubit - will template and rotate)
  constexpr double dp_cornerVertFudgeFactor = 60*CLHEP::um;
  constexpr double dp_cornerVertFudgeFactor2 = 1 * CLHEP::um; //Need this to ensure we don't get overlaps... (i.e. lengthen the base layer but not the flux line itself)
  constexpr double dp_cornerFluxLineBaseNbLayerDimX = 2105.000 * CLHEP::um; //was 2102.971
  constexpr double dp_cornerFluxLineBaseNbLayerDimY = 1948.084 * CLHEP::um + dp_cornerVertFudgeFactor + dp_cornerVertFudgeFactor2;
  constexpr double dp_cornerFluxLineBaseNbLayerDimZ = dp_groundPlaneDimZ;

  
  constexpr double cosOf45 = 0.7071067812; 
  constexpr double dp_cornerFluxLineCurveRadius = 100 * CLHEP::um;
  constexpr double dp_cornerFluxLineCurve1WrtPadPointX = cosOf45*dp_cornerFluxLineCurveRadius; //45 degree angle
  constexpr double dp_cornerFluxLineCurve1WrtPadPointY = cosOf45*dp_cornerFluxLineCurveRadius; //45 degree angle
  constexpr double dp_cornerFluxLineCurve1PadPointXWrtPadCenter = cosOf45*(0.5*dp_padEmptyPart1DimX + dp_padEmptyPart2TrdZ); //45 degree angle
  constexpr double dp_cornerFluxLineCurve1PadPointYWrtPadCenter = -1*cosOf45*(0.5*dp_padEmptyPart1DimX + dp_padEmptyPart2TrdZ); //45 degree angle
  constexpr double dp_cornerFluxLineCurve1WrtPadCenterX = dp_cornerFluxLineCurve1WrtPadPointX + dp_cornerFluxLineCurve1PadPointXWrtPadCenter;
  constexpr double dp_cornerFluxLineCurve1WrtPadCenterY = dp_cornerFluxLineCurve1WrtPadPointY + dp_cornerFluxLineCurve1PadPointYWrtPadCenter;
  constexpr double dp_cornerFluxLineHorizontalEmptyDimX = 1180.9309*CLHEP::um;
  constexpr double dp_cornerFluxLineHorizontalEmptyDimY = dp_fluxLineEmptyDimX;
  constexpr double dp_cornerFluxLineHorizontalEmptyDimZ = dp_groundPlaneDimZ;
  constexpr double dp_cornerFluxLineHorizontalConductorDimX = dp_cornerFluxLineHorizontalEmptyDimX;
  constexpr double dp_cornerFluxLineHorizontalConductorDimY = dp_fluxLineConductorDimX;
  constexpr double dp_cornerFluxLineHorizontalConductorDimZ = dp_groundPlaneDimZ;
  constexpr double dp_cornerFluxLineVerticalEmptyDimX = dp_fluxLineEmptyDimX;
  constexpr double dp_cornerFluxLineVerticalEmptyDimY = 1077.31*CLHEP::um + dp_cornerVertFudgeFactor; //20 is an arbitrary add to make length match center one
  constexpr double dp_cornerFluxLineVerticalEmptyDimZ = dp_groundPlaneDimZ;
  constexpr double dp_cornerFluxLineVerticalConductorDimX = dp_fluxLineConductorDimX;
  constexpr double dp_cornerFluxLineVerticalConductorDimY = dp_cornerFluxLineVerticalEmptyDimY;
  constexpr double dp_cornerFluxLineVerticalConductorDimZ = dp_groundPlaneDimZ;


  /*
  constexpr double dp_fluxLineEmptyDimX = 22 * CLHEP::um;
  constexpr double dp_fluxLineEmptyDimY = 1140.498 * CLHEP::um;
  constexpr double dp_fluxLineEmptyDimZ = dp_groundPlaneDimZ;
  constexpr double dp_fluxLineLineOffsetY = 0.5*dp_fluxLineBaseNbLayerDimY - dp_padEmptyPart1DimX - dp_padEmptyPart2TrdZ - 0.5*dp_fluxLineEmptyDimY;
  constexpr double dp_fluxLineConductorDimX = 10.15* CLHEP::um;
  constexpr double dp_fluxLineConductorDimY = dp_fluxLineEmptyDimY;
  constexpr double dp_fluxLineConductorDimZ = dp_groundPlaneDimZ;
  */


  
  //Overall flux line location parameters
  constexpr double dp_topCenterFluxLineOffsetX = 0;
  constexpr double dp_topCenterFluxLineOffsetY = 3812.378*CLHEP::um - 0.5*dp_fluxLineBaseNbLayerDimY;
  constexpr double dp_cornerFluxLinePadCornerDistanceFromWall = 51*CLHEP::um;
  constexpr double dp_topLeftFluxLineOffsetX = -0.5 * dp_groundPlaneDimX + 0.5 * dp_cornerFluxLineBaseNbLayerDimX + dp_cornerFluxLinePadCornerDistanceFromWall;
  constexpr double dp_topLeftFluxLineOffsetY = 0.5 * dp_groundPlaneDimY - 0.5 * dp_cornerFluxLineBaseNbLayerDimY - dp_cornerFluxLinePadCornerDistanceFromWall;





  //Corner flux line info
  
  
  
}







#endif
